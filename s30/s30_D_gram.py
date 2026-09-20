#!/usr/bin/env python
"""s30/s30_D_gram.py -- lane D: THE FIELD CLASS'S GRAM MATRIX, EFFECTIVE RANK, AND ORACLE
COMBINATION CEILING (the consequence of S30-L5).

S30-L5 established that the 21 native-free displacement fields of `s29/s29_D_fields.py` are not
noise: 11 of 21 have fold CIs excluding zero against the correct signed null and the best is
+7.7 sigma. Each is individually worth 0.0195 A because the transfer function is
RMSD_prod * sqrt(1 - rho^2), which is quadratic near zero. The question that follows, and which
nobody has asked: WHAT rho DOES THE BEST COMBINATION OF THE 21 REACH?

If the fields were mutually orthogonal, a combination would reach sqrt(sum rho_i^2) ~ 0.33-0.36,
and the bound's 3.00 A threshold is rho = 0.358. So the arithmetic lands near the threshold and
the answer must be MEASURED.

THE CONTROL THAT DECIDES IT, AND WHY IT IS NOT OPTIONAL. The fields live in the rigid-body-removed
space of dimension 3n - 6, which averages 32.9 over the 126 dev targets. Projecting the direction
to the native onto ANY 21-dimensional subspace of a 33-dimensional space captures, by dimension
counting alone, E[rho^2] ~ 21/33 = 0.64, i.e. rho ~ 0.80. An ORACLE combination ceiling quoted
without that control is almost entirely an artefact of how few residues these peptides have.
So every combination number here is reported beside RANDOM SUBSPACES OF THE SAME DIMENSION,
drawn with the same rigid-body removal (`s27/s28_A2_local.py`, A2's own reference machinery).

WHAT IS MEASURED
  1. the 21x21 Gram matrix of the unit field directions, per target and aggregated; its
     eigen-spectrum, stable rank (trace/lmax) and effective rank (exp of the spectral entropy)
  2. ORACLE rho* = ||P_span(D) u|| / ||u||, the best any weighting can do WITH THE NATIVE IN HAND
     -- beside a random-subspace control of matched dimension and of matched effective rank
  3. ORACLE GLOBAL weights (one w for all 126 targets, fitted with the native) -- the ceiling of a
     deployable rule, which per-target weights are not (S29's incidental-parameter result)
  4. the LEAVE-FOLD-OUT global weighting: fitted on 4 pinned folds, evaluated on the 5th, against
     the best single field chosen the same way. This is the only native-free number here.

Every rho in (2) and (3) is ORACLE and is labelled so in every line of output.

    python s30/s30_D_gram.py run [--slice i/n] [--limit N]
    python s30/s30_D_gram.py analyse

Regeneration is checked: the cosines recomputed here are asserted equal to the ones stored in
`s29/results/s29_D_fields_rows.jsonl` to 1e-9, so this file measures the SAME fields S30-L5 priced.
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
from s15 import seed as SD                 # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A2_local as A2         # noqa: E402
from s29 import s29_D_fields as F          # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(RESULTS, "s30_D_gram")
OUT = os.path.join(RESULTS, "s30_D_gram.json")
SALT = "s30Dgram"
TOL = 1e-8
N_CTRL = 64          # random-subspace control draws per target

#: the 11 fields whose fold-clustered CI excludes zero against the CORRECT signed null (S30-L5).
#: Fixed here in code, taken from s29/results/s29_D_fields.json, NOT re-selected on this file's
#: own numbers -- an equal-weight arm over fields chosen by the outcome would be leakage.
SIG11 = ("CHAN_DISTPOT", "MSET_250", "MSET_150", "CHAN_RG_LAW", "CHAN_CONTACT", "CONS_TRIM",
         "CHAN_LEG", "MSET_500", "CHAN_CONTACT_LL", "CHAN_ENV", "PROJ")

#: THE TWO BASES, AND THEY ARE NOT INTERCHANGEABLE. Every cosine here is measured on the CA POINT
#: CLOUD, where production is ~3.048 A. The bound's published thresholds are on the BUILT CHAIN,
#: where production is 3.2105 A (S29-L44 / REPORT_S29 section 0 item 2) -- which is why the record
#: says "rho = 0.358 is needed for 3.00 A": 3.2105*sqrt(1-0.358^2) = 3.00. Applying a point-cloud
#: production RMSD to that threshold silently answers a different question, so both are printed.
CHAIN_PROD = 3.2105


# ============================================================================ field vectors
def field_vectors(pdb):
    """The 21 native-free displacement fields as VECTORS, built exactly as `s29_D_fields.field_row`.

    Returns names, D (F, n, 3) rigid-removed displacements, u (n, 3) the ORACLE direction to the
    native (rigid-removed), C0 the production structure, and the per-field cosines for the
    regeneration check.
    """
    from s24 import d_harness as H
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(pdb)
    n, k = int(cand.n), int(cand.k)
    key = RP.rng_for(pdb, "tiekey").random(k)
    dis = np.asarray(ch["DIS"], float)
    order = np.lexsort((key, dis))
    top = np.sort(order[:F.PROD_M])
    C0, _ = H.readout_uniform(cand, top)
    C0 = np.asarray(C0, float)
    u = A2.remove_rigid(A2.oracle_direction(C0, cand.nat_ca), C0)            # ORACLE
    W = I.superpose_batch(np.asarray(cand.W, float), C0)
    names, D, cos = [], [], {}

    def add(name, C1):
        d = A2.remove_rigid(np.asarray(C1, float) - C0, C0)
        r = A2.rms(d)
        names.append(name); D.append(d)
        cos[name] = A2.cosine(d, u) if r > 1e-12 else float("nan")

    for m in F.M_LADDER:
        Cm, _ = H.readout_uniform(cand, np.sort(order[:m]))
        add(f"MSET_{m}", Cm)
    P = I.pairwise_rmsd(W[top])
    add("MEDOID", W[top][int(I.medoid(P))])
    pr = I.project(C0, cand.seq, cand.fold)
    add("PROJ", np.asarray(pr["ca"], float))
    for X in F.CHANNELS:
        if X not in ch:
            continue
        E = RP.combine(ch, ["DIS", X])
        oX = np.lexsort((key, np.asarray(E, float)))
        CX, _ = H.readout_uniform(cand, np.sort(oX[:F.PROD_M]))
        add(f"CHAN_{X}", CX)
    Pall = I.pairwise_rmsd(W)
    cons = Pall.sum(1) / max(k - 1, 1)
    oc = np.lexsort((key, cons))
    Cc, _ = H.readout_uniform(cand, np.sort(oc[:F.PROD_M]))
    add("CONS_TRIM", Cc)
    rg_pool = float(np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1)).mean())
    cen = C0.mean(0)
    rg0 = float(np.sqrt(((C0 - cen) ** 2).sum(1).mean()))
    add("EXPAND", cen + (C0 - cen) * (rg_pool / max(rg0, 1e-12)))
    return dict(pdb=pdb, n=n, fold=int(cand.fold), fail18=bool(pdb in I.FAIL18), names=names,
                D=np.stack(D), u=u, C0=C0, cos=cos,
                rmsd_prod=float(I.ca_rmsd(C0, cand.nat_ca)))


_S29 = None


def s29_cos(pdb):
    global _S29
    if _S29 is None:
        _S29 = {}
        p = os.path.join(ROOT, "s29", "results", "s29_D_fields_rows.jsonl")
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line); _S29[r["pdb"]] = r["cos"]
    return _S29.get(pdb)


def build(pdbs):
    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    for q, pdb in enumerate(pdbs):
        f = os.path.join(CACHE, pdb + ".npz")
        if os.path.exists(f):
            continue
        t1 = time.time()
        r = field_vectors(pdb)
        ref = s29_cos(pdb)
        if ref:                                         # implementation-drift check
            dev = max(abs(r["cos"][k] - ref[k]) for k in r["names"] if k in ref and np.isfinite(ref[k]))
            assert dev < 1e-9, f"{pdb}: regenerated cosines differ from S29's by {dev:.2e}"
        tmp = f + f".tmp{os.getpid()}.npz"
        np.savez_compressed(tmp, D=r["D"], u=r["u"], C0=r["C0"], names=np.array(r["names"]),
                            n=np.array(r["n"]), fold=np.array(r["fold"]), fail18=np.array(int(r["fail18"])),
                            rmsd_prod=np.array(r["rmsd_prod"]),
                            cos=np.array([r["cos"][k] for k in r["names"]]))
        from s27 import s28_A_amp as A
        A.replace_retry(tmp, f)
        print(f"  [{q+1}/{len(pdbs)}] {pdb} n={r['n']} F={len(r['names'])} cosines match S29 "
              f"({time.time()-t1:.1f}s, {(time.time()-t0)/60:.1f} min)", flush=True)
    print("cache:", CACHE)


def load(pdb):
    f = os.path.join(CACHE, pdb + ".npz")
    if not os.path.exists(f):
        return None
    with np.load(f, allow_pickle=False) as z:
        return dict(pdb=pdb, D=np.array(z["D"]), u=np.array(z["u"]), C0=np.array(z["C0"]),
                    names=[str(x) for x in z["names"]], n=int(z["n"]), fold=int(z["fold"]),
                    fail18=bool(int(z["fail18"])), rmsd_prod=float(z["rmsd_prod"]),
                    cos=np.array(z["cos"]))


# ============================================================================ geometry
def _flat_unit(D):
    """(F, 3n) unit-norm rows; rows with zero norm dropped, their indices returned."""
    X = D.reshape(len(D), -1)
    nrm = np.linalg.norm(X, axis=1)
    keep = nrm > 1e-12
    return X[keep] / nrm[keep, None], keep


def proj_rho(X, u):
    """||P_rowspace(X) u|| / ||u||: the best cosine ANY weighting of the rows of X achieves."""
    uu = u.ravel()
    nu = np.linalg.norm(uu)
    if nu < 1e-15 or len(X) == 0:
        return float("nan")
    Q, _ = np.linalg.qr(X.T)                     # orthonormal basis of the row space
    return float(np.linalg.norm(Q.T @ uu) / nu)


def rand_subspace_rho(C0, u, dim, n_draw, pdb):
    """The same projection onto `dim` RANDOM directions, rigid body removed exactly as A2 does."""
    rng = SD.stable_rng(pdb, "s30Dgram_ctrl", salt=SALT)
    out = []
    for _ in range(n_draw):
        R = np.stack([A2.remove_rigid(rng.normal(size=C0.shape), C0) for _ in range(dim)])
        X, _ = _flat_unit(R)
        out.append(proj_rho(X, u))
    return np.array(out, float)


def spectrum(X):
    """Gram of unit rows, its eigenvalues, stable rank and effective rank."""
    G = X @ X.T
    w = np.linalg.eigvalsh(G)[::-1]
    w = np.clip(w, 0, None)
    tr = float(w.sum())
    lmax = float(w[0]) if len(w) else float("nan")
    p = w / tr if tr > 0 else w
    ent = float(-np.sum(p[p > 0] * np.log(p[p > 0])))
    return dict(eig=[float(x) for x in w], trace=tr, lmax=lmax,
                stable_rank=float(tr / lmax) if lmax > 0 else float("nan"),
                effective_rank=float(np.exp(ent)),
                numerical_rank=int((w > TOL * lmax).sum()),
                G=G)


# ============================================================================ analysis
def analyse(pdbs):
    rows = []
    for pdb in pdbs:
        r = load(pdb)
        if r is None:
            continue
        X, keep = _flat_unit(r["D"])
        names = [nm for nm, k in zip(r["names"], keep) if k]
        sp = spectrum(X)
        amb = 3 * r["n"] - 6
        rho_or = proj_rho(X, r["u"])
        rr = int(sp["numerical_rank"])
        er = int(round(sp["effective_rank"]))
        ctrl_full = rand_subspace_rho(r["C0"], r["u"], min(len(X), amb), N_CTRL, r["pdb"])
        ctrl_eff = rand_subspace_rho(r["C0"], r["u"], max(1, min(er, amb)), N_CTRL, r["pdb"])
        rows.append(dict(pdb=r["pdb"], n=r["n"], fold=r["fold"], fail18=r["fail18"], amb=amb,
                         rmsd_prod=r["rmsd_prod"], names=names, F=len(X),
                         stable_rank=sp["stable_rank"], effective_rank=sp["effective_rank"],
                         numerical_rank=rr, eig=sp["eig"], G=sp["G"].tolist(),
                         rho_oracle=rho_or, best_single=float(np.nanmax(np.abs(r["cos"]))),
                         ctrl_full_mean=float(np.nanmean(ctrl_full)), ctrl_full_p95=float(np.nanpercentile(ctrl_full, 95)),
                         ctrl_eff_mean=float(np.nanmean(ctrl_eff)), ctrl_eff_p95=float(np.nanpercentile(ctrl_eff, 95)),
                         cos=[float(x) for x in r["cos"]]))
    if not rows:
        raise SystemExit("no cached targets; run `python s30/s30_D_gram.py run` first")
    names0 = rows[0]["names"]
    ok = [r for r in rows if r["names"] == names0]
    pdbl = [r["pdb"] for r in ok]
    folds = ST.pinned_folds(pdbl)
    fail = np.array([r["fail18"] for r in ok], bool)
    rp = float(np.mean([r["rmsd_prod"] for r in ok]))

    G = np.mean([np.array(r["G"]) for r in ok], 0)
    Gsp = np.linalg.eigvalsh(G)[::-1]
    out = dict(n=len(ok), n_all=len(rows), fields=names0, rmsd_prod_mean=rp,
               ambient_mean=float(np.mean([r["amb"] for r in ok])),
               oracle=True, salt=SALT, n_ctrl=N_CTRL)
    out["gram_mean"] = dict(matrix=G.tolist(), eig=[float(x) for x in Gsp],
                            stable_rank=float(Gsp.sum() / Gsp[0]),
                            effective_rank=float(np.exp(-np.sum((Gsp / Gsp.sum()) * np.log(np.clip(Gsp / Gsp.sum(), 1e-300, None))))),
                            top3_share=float(Gsp[:3].sum() / Gsp.sum()),
                            mean_abs_offdiag=float(np.abs(G - np.diag(np.diag(G)))[np.triu_indices(len(G), 1)].mean()))
    out["per_target_rank"] = dict(
        stable=ST.compare(np.array([r["stable_rank"] for r in ok]), np.zeros(len(ok)), folds, seed_parts=(SALT,), label="stable rank")["effect"],
        effective=float(np.mean([r["effective_rank"] for r in ok])),
        numerical=float(np.mean([r["numerical_rank"] for r in ok])),
        numerical_min=int(min(r["numerical_rank"] for r in ok)),
        numerical_max=int(max(r["numerical_rank"] for r in ok)))

    # --- (2) the ORACLE per-target combination, against matched random subspaces
    ro = np.array([r["rho_oracle"] for r in ok]); cf = np.array([r["ctrl_full_mean"] for r in ok])
    ce = np.array([r["ctrl_eff_mean"] for r in ok]); bs = np.array([r["best_single"] for r in ok])
    c1 = ST.compare(ro, cf, folds, names=pdbl, seed_parts=(SALT,), label="ORACLE per-target rho* - matched random 21-dim subspace")
    out["oracle_per_target"] = dict(
        rho=float(ro.mean()), rho_median=float(np.median(ro)),
        ctrl_full_dim_rho=float(cf.mean()), ctrl_eff_dim_rho=float(ce.mean()),
        best_single_abs_cos=float(bs.mean()),
        implied_rmsd=rp * float(np.sqrt(max(1 - ro.mean() ** 2, 0))),
        implied_rmsd_ctrl=rp * float(np.sqrt(max(1 - cf.mean() ** 2, 0))),
        excess_over_control={k: v for k, v in c1.items() if k != "concentration"}, fmt=ST.fmt(c1))

    # --- (3) ORACLE GLOBAL weights (one w, all targets, native in hand) and (4) leave-fold-out
    Xs, us, C0s = [], [], []
    for r in ok:
        z = load(r["pdb"])
        X, keep = _flat_unit(z["D"])
        uu = z["u"].ravel(); uu = uu / np.linalg.norm(uu)
        Xs.append(X); us.append(uu); C0s.append(z["C0"])

    # --- (2b) THE RANK CURVE. The unconstrained ORACLE projection uses all 21 numerically
    #     independent directions, which in a ~33-dimensional space is mostly dimension counting
    #     and needs unbounded weights on near-degenerate modes. The informative object is the
    #     projection onto the top r PRINCIPAL directions of the field set, r = 1..F, each read
    #     beside a RANDOM r-dimensional subspace of the same space.
    F = len(names0)
    rk_real = np.full((len(ok), F), np.nan); rk_ctrl = np.full((len(ok), F), np.nan)
    for t, r in enumerate(ok):
        U, S, Vt = np.linalg.svd(Xs[t], full_matrices=False)
        uu = us[t]
        for rr in range(1, F + 1):
            B = Vt[:rr]
            rk_real[t, rr - 1] = float(np.linalg.norm(B @ uu))
        rng = SD.stable_rng(r["pdb"], "s30Dgram_rankctrl", salt=SALT)
        R = np.stack([A2.remove_rigid(rng.normal(size=C0s[t].shape), C0s[t]).ravel() for _ in range(F)])
        Qr, _ = np.linalg.qr(R.T)
        for rr in range(1, F + 1):
            rk_ctrl[t, rr - 1] = float(np.linalg.norm(Qr[:, :rr].T @ uu))
    out["rank_curve"] = dict(
        r=list(range(1, F + 1)),
        rho_real=[float(np.nanmean(rk_real[:, j])) for j in range(F)],
        rho_ctrl=[float(np.nanmean(rk_ctrl[:, j])) for j in range(F)],
        excess=[float(np.nanmean(rk_real[:, j] - rk_ctrl[:, j])) for j in range(F)],
        excess_over_mde=[float(ST.compare(rk_real[:, j], rk_ctrl[:, j], folds, seed_parts=(SALT,),
                                          label=f"rank {j+1}")["effect_over_mde"]) for j in range(F)],
        ci95_fold=[ST.compare(rk_real[:, j], rk_ctrl[:, j], folds, seed_parts=(SALT,),
                              label=f"rank {j+1}")["ci95_fold"] for j in range(F)],
        implied_rmsd_real=[rp * float(np.sqrt(max(1 - np.nanmean(rk_real[:, j]) ** 2, 0))) for j in range(F)],
        implied_rmsd_ctrl=[rp * float(np.sqrt(max(1 - np.nanmean(rk_ctrl[:, j]) ** 2, 0))) for j in range(F)],
        note="rho at rank r for the REAL field set's top-r principal directions vs a RANDOM r-dimensional "
             "subspace of the same rigid-body-removed space. Both are ORACLE (u is the direction to the "
             "native). The excess is the only part of the ceiling the fields themselves earn.")

    LAMS = (1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e5)

    def fit_global(idx, lam):
        A = np.zeros((len(names0), len(names0))); b = np.zeros(len(names0))
        for t in idx:
            A += Xs[t] @ Xs[t].T; b += Xs[t] @ us[t]
        return np.linalg.solve(A + lam * np.trace(A) / len(A) * np.eye(len(A)), b)

    def cos_of(w, t):
        v = Xs[t].T @ w
        nv = np.linalg.norm(v)
        return float(v @ us[t] / nv) if nv > 1e-15 else 0.0

    #: cos_t(w) = (w.a_t) / sqrt(w' G_t w) -- exact, so the ORACLE-optimal global weighting can be
    #: found by MAXIMISING THE MEAN COSINE rather than by least squares. The LS solution is NOT the
    #: optimum of this objective: on the partial run it returned rho BELOW the best single field,
    #: which is impossible for a true maximum (w = e_j reproduces field j exactly). Caught by that
    #: internal consistency check, which is why the check is worth running.
    Amat = np.stack([Xs[t] @ us[t] for t in range(len(ok))])                     # (T, F)
    Gmat = np.stack([Xs[t] @ Xs[t].T for t in range(len(ok))])                   # (T, F, F)

    def cos_vec(w, idx):
        a = Amat[idx] @ w
        s = np.sqrt(np.maximum(np.einsum("tij,i,j->t", Gmat[idx], w, w), 1e-300))
        return a / s

    def neg_mean_cos(w, idx):
        a = Amat[idx] @ w
        Gw = np.einsum("tij,j->ti", Gmat[idx], w)
        s = np.sqrt(np.maximum((Gw * w).sum(1), 1e-300))
        c = a / s
        g = (Amat[idx] / s[:, None]) - (a / s ** 3)[:, None] * Gw
        return -float(c.mean()), -g.mean(0)

    def max_cos(idx, extra_starts=()):
        """The ORACLE-optimal GLOBAL weighting: argmax_w mean_t cos(X_t' w, u_t), multi-start.

        Starts: every single-field basis vector (so the optimum can never be worse than the best
        single field -- the consistency floor), equal weights, and any extra start supplied.
        """
        from scipy.optimize import minimize
        idx = np.asarray(idx, int)
        starts = [np.eye(len(names0))[j] for j in range(len(names0))]
        starts.append(np.ones(len(names0)))
        starts.extend(extra_starts)
        best_w, best_v = None, -np.inf
        for w0 in starts:
            r = minimize(neg_mean_cos, np.asarray(w0, float), args=(idx,), jac=True,
                         method="L-BFGS-B", options=dict(maxiter=500))
            v = -float(r.fun)
            if v > best_v:
                best_v, best_w = v, r.x
        return best_w, best_v

    def pick_lam(idx):
        """The ridge chosen by NESTED leave-one-fold-out INSIDE the training folds only.
        Never on the held-out fold: the collinearity here is severe (stable rank ~2), so an
        untuned ridge would be a strawman and a test-tuned one would be leakage."""
        sub = np.array([folds[t] for t in idx])
        best, bl = -np.inf, LAMS[0]
        for lam in LAMS:
            v = []
            for q in sorted(set(sub.tolist())):
                tr2 = [t for t, f in zip(idx, sub) if f != q]
                te2 = [t for t, f in zip(idx, sub) if f == q]
                if not tr2 or not te2:
                    continue
                w = fit_global(tr2, lam)
                v += [cos_of(w, t) for t in te2]
            m = float(np.mean(v)) if v else -np.inf
            if m > best:
                best, bl = m, lam
        return bl, best

    C = np.array([r["cos"] for r in ok], float)
    lam_all, _ = pick_lam(list(range(len(ok))))
    wall_ls = fit_global(range(len(ok)), lam_all)
    wall, _ = max_cos(range(len(ok)), extra_starts=(wall_ls,))            # the TRUE oracle optimum
    cos_global_oracle = np.array([cos_of(wall, t) for t in range(len(ok))])
    cos_global_ls = np.array([cos_of(wall_ls, t) for t in range(len(ok))])
    lfo = np.zeros(len(ok)); lfo_single = np.zeros(len(ok)); lfo_eq = np.zeros(len(ok)); lams_used = {}
    sel_used = {}
    for q in sorted(set(folds.tolist())):
        tr = [t for t in range(len(ok)) if folds[t] != q]
        te = [t for t in range(len(ok)) if folds[t] == q]
        lam, _ = pick_lam(tr)
        lams_used[int(q)] = lam
        w_ls = fit_global(tr, lam)
        w, _ = max_cos(tr, extra_starts=(w_ls,))       # fitted on TRAIN ONLY, same objective
        for t in te:
            lfo[t] = cos_of(w, t)
        j = int(np.nanargmax(np.nanmean(C[tr], 0)))          # best single field on the TRAIN folds
        for t in te:
            lfo_single[t] = C[t, j]
        #: equal weights over the fields SELECTED INSIDE THE TRAINING FOLDS -- the honest version
        #: of EQ11, whose 11 were chosen on all 126 targets and is therefore mildly leaked.
        mtr = np.nanmean(C[tr], 0)
        sd = np.nanstd(C[tr], 0, ddof=1) / max(np.sqrt(len(tr)), 1)
        pick = np.where(mtr > 1.96 * sd)[0]
        if len(pick) == 0:
            pick = np.array([j])
        sel_used[int(q)] = [names0[i] for i in pick]
        wq = np.zeros(len(names0)); wq[pick] = 1.0
        for t in te:
            lfo_eq[t] = cos_of(wq, t)
    #: two ZERO-PARAMETER arms: nothing is fitted, so nothing can be overfitted. The honest
    #: simple control for a 21-parameter rule (`does it survive a simpler control?`).
    eq_all = np.array([cos_of(np.ones(len(names0)), t) for t in range(len(ok))])
    sig = [i for i, nm in enumerate(names0) if nm in SIG11]
    wsig = np.zeros(len(names0)); wsig[sig] = 1.0
    eq_sig = np.array([cos_of(wsig, t) for t in range(len(ok))])
    c2 = ST.compare(lfo, lfo_single, folds, names=pdbl, seed_parts=(SALT,), label="LFO global weighting - LFO best single field")
    c3 = ST.compare(lfo, np.zeros(len(ok)), folds, names=pdbl, seed_parts=(SALT,), label="LFO global weighting vs 0")
    c4 = ST.compare(eq_sig, lfo_single, folds, names=pdbl, seed_parts=(SALT,), label="equal-weight 11 significant fields - LFO best single field")
    c5 = ST.compare(eq_sig, np.zeros(len(ok)), folds, names=pdbl, seed_parts=(SALT,), label="equal-weight 11 vs 0")
    c6 = ST.compare(lfo_eq, lfo_single, folds, names=pdbl, seed_parts=(SALT,), label="LFO-SELECTED equal weights - LFO best single field")
    c7 = ST.compare(lfo_eq, np.zeros(len(ok)), folds, names=pdbl, seed_parts=(SALT,), label="LFO-SELECTED equal weights vs 0")
    out["global_weights"] = dict(
        oracle_global_rho=float(cos_global_oracle.mean()), oracle_global_lambda=lam_all,
        oracle_global_implied_rmsd=rp * float(np.sqrt(max(1 - cos_global_oracle.mean() ** 2, 0))),
        oracle_global_rho_ls=float(cos_global_ls.mean()),
        oracle_global_best_single_floor=float(np.nanmax(np.nanmean(C, 0))),
        lfo_rho=float(lfo.mean()), lfo_median=float(np.median(lfo)), lfo_lambdas=lams_used,
        lfo_implied_rmsd=rp * float(np.sqrt(max(1 - lfo.mean() ** 2, 0))),
        lfo_best_single_rho=float(lfo_single.mean()),
        equal_all_rho=float(eq_all.mean()), equal_sig11_rho=float(eq_sig.mean()),
        equal_sig11_implied_rmsd=rp * float(np.sqrt(max(1 - eq_sig.mean() ** 2, 0))),
        lfo_eq_rho=float(lfo_eq.mean()), lfo_eq_median=float(np.median(lfo_eq)),
        lfo_eq_implied_rmsd=rp * float(np.sqrt(max(1 - lfo_eq.mean() ** 2, 0))),
        lfo_eq_selected=sel_used,
        lfo_eq_vs_zero={k: v for k, v in c7.items() if k != "concentration"},
        lfo_eq_vs_best_single={k: v for k, v in c6.items() if k != "concentration"},
        lfo_eq_concentration=c7.get("concentration"), fmt_lfo_eq=ST.fmt(c6),
        lfo_vs_zero={k: v for k, v in c3.items() if k != "concentration"},
        lfo_vs_best_single={k: v for k, v in c2.items() if k != "concentration"},
        eq11_vs_best_single={k: v for k, v in c4.items() if k != "concentration"},
        eq11_vs_zero={k: v for k, v in c5.items() if k != "concentration"},
        fmt_vs_zero=ST.fmt(c3), fmt_vs_single=ST.fmt(c2), fmt_eq11=ST.fmt(c4),
        weights={nm: float(v) for nm, v in zip(names0, wall)})
    out["fail18_split"] = dict(
        rho_oracle=dict(fail=float(ro[fail].mean()), other=float(ro[~fail].mean())),
        lfo=dict(fail=float(lfo[fail].mean()), other=float(lfo[~fail].mean())))
    def chain(r):
        return CHAIN_PROD * float(np.sqrt(max(1 - r ** 2, 0)))
    out["thresholds"] = dict(
        cloud_prod=rp, chain_prod=CHAIN_PROD,
        rho_for_3A_cloud=float(np.sqrt(max(1 - (3.00 / rp) ** 2, 0))),
        rho_for_3A_chain=float(np.sqrt(max(1 - (3.00 / CHAIN_PROD) ** 2, 0))),
        rho_for_2_5A_chain=float(np.sqrt(max(1 - (2.50 / CHAIN_PROD) ** 2, 0))),
        note="rho needed so that rmsd_prod*sqrt(1-rho^2) reaches the stated RMSD. The record's "
             "'0.358 for 3.00 A' is the CHAIN row; every cosine measured in this file is a POINT-CLOUD "
             "cosine, so the chain column below is the comparable one.")
    out["chain_basis"] = dict(
        oracle_per_target=chain(out["oracle_per_target"]["rho"]), oracle_global=chain(out["global_weights"]["oracle_global_rho"]),
        lfo_global=chain(out["global_weights"]["lfo_rho"]), lfo_eq=chain(out["global_weights"]["lfo_eq_rho"]),
        eq_sig11=chain(out["global_weights"]["equal_sig11_rho"]), best_single=chain(out["global_weights"]["oracle_global_best_single_floor"]),
        ctrl_random_subspace=chain(out["oracle_per_target"]["ctrl_full_dim_rho"]), production=CHAIN_PROD,
        note="the SAME cosines carried through the bound with the BUILT-CHAIN production RMSD. "
             "These are the numbers comparable to the record's 3.2105 / 3.00 / 2.50 A.")
    out["text"] = render(out)
    os.makedirs(RESULTS, exist_ok=True)
    ST.save_atomic(OUT, out, module_file=__file__)
    print(out["text"])
    print("wrote", OUT)
    return out


def render(o):
    g = o["gram_mean"]; pr = o["per_target_rank"]; op = o["oracle_per_target"]; gw = o["global_weights"]
    L = [f"FIELD-CLASS GRAM AND COMBINATION CEILING (ORACLE throughout), n={o['n']} targets, "
         f"{len(o['fields'])} fields, mean ambient dim (3n-6) {o['ambient_mean']:.1f}, production {o['rmsd_prod_mean']:.4f} A"]
    th = o["thresholds"]
    L.append(f"  BASES: cosines are POINT-CLOUD (production {th['cloud_prod']:.4f} A). The record's thresholds are BUILT CHAIN "
             f"(production {th['chain_prod']:.4f} A): rho for 3.00 A = {th['rho_for_3A_chain']:.4f}, for 2.50 A = {th['rho_for_2_5A_chain']:.4f}.")
    L.append("  (1) GRAM of the unit field directions, averaged over targets:")
    L.append(f"      mean |off-diagonal| {g['mean_abs_offdiag']:.4f}  eigenvalues {['%.2f' % x for x in g['eig'][:8]]} ...")
    L.append(f"      stable rank (trace/lmax) {g['stable_rank']:.3f}   effective rank (exp entropy) {g['effective_rank']:.3f}   "
             f"top-3 eigenvalue share {100*g['top3_share']:.1f}%")
    L.append(f"      per-target: stable {pr['stable']:.3f}  effective {pr['effective']:.3f}  numerical {pr['numerical']:.2f} "
             f"(min {pr['numerical_min']}, max {pr['numerical_max']}) out of {len(o['fields'])} fields in {o['ambient_mean']:.0f} dimensions")
    L.append("  (2) ORACLE per-target combination (best weighting WITH THE NATIVE IN HAND), vs matched random subspaces:")
    L.append(f"      ORACLE rho* {op['rho']:.4f} (median {op['rho_median']:.4f})  -> implied {op['implied_rmsd']:.4f} A")
    L.append(f"      matched RANDOM subspace, same dimension: {op['ctrl_full_dim_rho']:.4f} -> implied {op['implied_rmsd_ctrl']:.4f} A")
    L.append(f"      matched RANDOM subspace, effective rank:  {op['ctrl_eff_dim_rho']:.4f}")
    L.append(f"      best SINGLE field |cos| (oracle-chosen per target): {op['best_single_abs_cos']:.4f}")
    e = op["excess_over_control"]
    L.append(f"      EXCESS over the matched random subspace: {e['effect']:+.4f}  SE {e['se']:.4f}  {e['effect_over_mde']:+.2f}x MDE  "
             f"fold CI [{e['ci95_fold'][0]:+.3f}, {e['ci95_fold'][1]:+.3f}]  folds same sign {e['folds_same_sign']}/{e['n_folds']}")
    rc = o.get("rank_curve")
    if rc:
        L.append("  (2b) THE RANK CURVE -- ORACLE rho at rank r (top-r principal directions of the field set) "
                 "vs a RANDOM r-dimensional subspace of the same space:")
        L.append("        r   rho_real  rho_rand    excess   xMDE   fold CI            implied A (real / rand)")
        for j, r in enumerate(rc["r"]):
            ci = rc["ci95_fold"][j]
            L.append(f"      {r:3d}   {rc['rho_real'][j]:8.4f}  {rc['rho_ctrl'][j]:8.4f}  {rc['excess'][j]:+8.4f}  "
                     f"{rc['excess_over_mde'][j]:+5.2f}  [{ci[0]:+.3f},{ci[1]:+.3f}]   {rc['implied_rmsd_real'][j]:.4f} / {rc['implied_rmsd_ctrl'][j]:.4f}")
    L.append("  (3) ORACLE GLOBAL weighting (ONE w for all targets, fitted with the native by MAXIMISING the mean cosine):")
    L.append(f"      rho {gw['oracle_global_rho']:.4f} -> implied {gw['oracle_global_implied_rmsd']:.4f} A   "
             f"| consistency floor (best single field's own mean cos) {gw['oracle_global_best_single_floor']:+.4f}  "
             f"| least-squares fit for comparison {gw['oracle_global_rho_ls']:.4f} (ridge {gw['oracle_global_lambda']:g})")
    L.append("  (4) NATIVE-FREE ARMS -- the only numbers here that may be read as a capability:")
    L.append(f"      LFO global weighting (ridge by NESTED CV inside train, lambdas {gw['lfo_lambdas']}):")
    L.append(f"        rho {gw['lfo_rho']:.4f} (median {gw['lfo_median']:.4f}) -> implied {gw['lfo_implied_rmsd']:.4f} A")
    L.append(f"      LFO best SINGLE field (the simpler control): rho {gw['lfo_best_single_rho']:.4f}")
    L.append(f"      ZERO-PARAMETER equal weights over the 11 S30-L5-significant fields: rho {gw['equal_sig11_rho']:.4f} "
             f"-> implied {gw['equal_sig11_implied_rmsd']:.4f} A   | over all 21: {gw['equal_all_rho']:.4f}")
    z = gw["lfo_vs_zero"]; s = gw["lfo_vs_best_single"]; e1 = gw["eq11_vs_zero"]; e2 = gw["eq11_vs_best_single"]
    L.append(f"      LFO vs zero:             {z['effect']:+.4f}  {z['effect_over_mde']:+.2f}x MDE  fold CI [{z['ci95_fold'][0]:+.3f}, {z['ci95_fold'][1]:+.3f}]  folds {z['folds_same_sign']}/{z['n_folds']}")
    L.append(f"      LFO vs best single fld:  {s['effect']:+.4f}  {s['effect_over_mde']:+.2f}x MDE  fold CI [{s['ci95_fold'][0]:+.3f}, {s['ci95_fold'][1]:+.3f}]  folds {s['folds_same_sign']}/{s['n_folds']}")
    L.append(f"      EQ11 vs zero:            {e1['effect']:+.4f}  {e1['effect_over_mde']:+.2f}x MDE  fold CI [{e1['ci95_fold'][0]:+.3f}, {e1['ci95_fold'][1]:+.3f}]  folds {e1['folds_same_sign']}/{e1['n_folds']}")
    L.append(f"      EQ11 vs best single fld: {e2['effect']:+.4f}  {e2['effect_over_mde']:+.2f}x MDE  fold CI [{e2['ci95_fold'][0]:+.3f}, {e2['ci95_fold'][1]:+.3f}]  folds {e2['folds_same_sign']}/{e2['n_folds']}")
    f1 = gw["lfo_eq_vs_zero"]; f2 = gw["lfo_eq_vs_best_single"]
    L.append(f"      LFO-SELECTED equal weights (fields chosen INSIDE train, no leakage): rho {gw['lfo_eq_rho']:.4f} "
             f"(median {gw['lfo_eq_median']:.4f}) -> implied {gw['lfo_eq_implied_rmsd']:.4f} A")
    L.append(f"        vs zero:            {f1['effect']:+.4f}  {f1['effect_over_mde']:+.2f}x MDE  fold CI [{f1['ci95_fold'][0]:+.3f}, {f1['ci95_fold'][1]:+.3f}]  folds {f1['folds_same_sign']}/{f1['n_folds']}")
    L.append(f"        vs best single fld: {f2['effect']:+.4f}  {f2['effect_over_mde']:+.2f}x MDE  fold CI [{f2['ci95_fold'][0]:+.3f}, {f2['ci95_fold'][1]:+.3f}]  folds {f2['folds_same_sign']}/{f2['n_folds']}")
    L.append(f"        fields selected per fold: {gw['lfo_eq_selected']}")
    cc = gw.get("lfo_eq_concentration")
    if cc:
        L.append(f"        concentration vs the UNIFORM-EFFECT null: median {cc['median']:+.4f} drop-top10 {cc['drop_top10_mean']:+.4f} "
                 f"at the {100*cc['pctile_in_null']:.0f}th pct of the null, flag {cc['flag']}")
    L.append(f"  (5) FAIL18 vs 108: ORACLE rho* {o['fail18_split']['rho_oracle']['fail']:.4f} / {o['fail18_split']['rho_oracle']['other']:.4f}   "
             f"LFO {o['fail18_split']['lfo']['fail']:.4f} / {o['fail18_split']['lfo']['other']:.4f}")
    cb = o["chain_basis"]
    L.append(f"  (6) THE SAME COSINES ON THE BUILT CHAIN (production {cb['production']:.4f} A) -- the comparable column:")
    L.append(f"      ORACLE per-target {cb['oracle_per_target']:.4f}  (matched random subspace {cb['ctrl_random_subspace']:.4f})  "
             f"| ORACLE GLOBAL {cb['oracle_global']:.4f}  | best single field {cb['best_single']:.4f}")
    L.append(f"      native-free: LFO global {cb['lfo_global']:.4f}   LFO-selected equal weights {cb['lfo_eq']:.4f}   "
             f"EQ11 (leaked selection) {cb['eq_sig11']:.4f}")
    L.append("  READ (2) ONLY BESIDE ITS CONTROL: a 21-dim subspace of a ~33-dim space captures most of any direction "
             "by dimension counting. The ORACLE ceiling is meaningful only as its EXCESS over the matched random subspace.")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "analyse"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--slice", default=None)
    a = ap.parse_args()
    from s25 import phys_lib as P
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    if a.slice:
        i, n = (int(x) for x in a.slice.split("/"))
        pdbs = [p for q, p in enumerate(pdbs) if q % n == i]
    if a.mode == "run":
        build(pdbs)
    else:
        analyse(pdbs)


if __name__ == "__main__":
    main()
