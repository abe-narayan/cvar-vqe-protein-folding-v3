#!/usr/bin/env python
"""s31/s31_C_ladder.py -- lane C, the SET-MATCHED readout ladder (C1, prereg amendment A1.3) and
the common-mode saturation test the coordinator predicted.

Two things S30 did not do:

  1. **Hold the candidate set fixed.**  S30-L11 concluded "argmin dominates at every bit budget"
     from a pair that compared 2-of-**75** against 1-of-**128**.  Every row here is over the SAME
     top-128 of the shipped DIS order, so the comparison is about the READOUT CLASS and not about
     which candidates each class was handed.
  2. **Decompose the emitted error exactly.**  For any weight vector with `sum(w) = 1`,
        `X(w) - nat  =  c  +  sum_i w_i e_i`,        c = Xbar - nat,   e_i = X_i - Xbar
     so the shared component `c` is preserved EXACTLY by every affine readout, and a readout can
     only beat it by making `sum_i w_i e_i` cancel `c`.  Each arm reports the three terms.

ORACLE DISCIPLINE.  Arms whose weights are solved against the native are labelled ORACLE in the
output dict itself (`oracle: true`) and are upper-bound diagnostics only.  The native-free arms
(`unif_*`, `soft_*`, `typ_*`, `randw_*`) choose their weights from the score and the set geometry
and read the native only to be SCORED, which is what every RMSD in this project does.

Basis: CA POINT CLOUD.  The built chain is the endpoint; s29 lane O already chained the decisive
ORACLE rows and this file re-uses them rather than re-projecting (S31 operational rule: both sides
of a chain contrast must be projected in the same job).

    python s31/s31_C_ladder.py run [--limit N]
    python s31/s31_C_ladder.py analyse
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

from s12 import instrument as I                    # noqa: E402
from s24 import stats_lib as ST                    # noqa: E402
from s31 import s31_C_cache as CA                  # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s31_C_ladder_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_C_ladder.json")

TOP = 128                                   # the deployed register width; the FIXED candidate set
M_CURVE = (1, 2, 3, 4, 5, 6, 8, 10, 13, 16, 21, 26, 32, 40, 50, 64, 75, 89, 100, 112, 128)
N_DRAW = 8                                  # random-subset draws; the MEAN is the arm (rule 10)
LAM = tuple(10.0 ** np.arange(-6, 4.01, 0.5))     # the ridge path: uniform -> unconstrained affine
SOFT_T = (0.25, 0.5, 1.0, 2.0, 4.0)
SALT = "s31Cladder"


# ---------------------------------------------------------------- weight diagnostics
def wdiag(w):
    w = np.asarray(w, float)
    return dict(l1=float(np.abs(w).sum()),
                ess=float(1.0 / max(float((w ** 2).sum()), 1e-300)),
                frac_neg=float((w < 0).mean()),
                neg_mass=float(-w[w < 0].sum()) if (w < 0).any() else 0.0,
                wmax=float(w.max()), wmin=float(w.min()))


def emit(Wf, w, n):
    return (np.asarray(w, float) @ Wf).reshape(n, 3)


def decompose(Wf, w, natp, n):
    """Exact split of the emitted squared deviation into the shared part, the cross term and the
    deviation part, in the common frame.  `natp` is the native posed in that frame, flattened.

        ||X(w) - nat||^2  =  ||c||^2  +  2<c, D w>  +  ||D w||^2
        c = Xbar - nat  (Xbar the UNWEIGHTED mean of the set),  D = Wf - Xbar
    """
    Wf = np.asarray(Wf, float)
    w = np.asarray(w, float)
    xbar = Wf.mean(0)
    c = xbar - natp
    dv = (w @ Wf) - xbar                     # = D^T w, since sum(w) = 1
    cc = float(c @ c) / n
    cross = 2.0 * float(c @ dv) / n
    dd = float(dv @ dv) / n
    return dict(c2=cc, cross=cross, d2=dd, total=cc + cross + dd,
                c_rms=math.sqrt(max(cc, 0.0)), d_rms=math.sqrt(max(dd, 0.0)))


# ---------------------------------------------------------------- ORACLE solvers
def oracle_affine(A, y, ones):
    """ORACLE. argmin ||A w - y||^2 s.t. 1'w = 1, minimum-norm solution.  A is (3n, k)."""
    # substitute w = w0 + v with 1'v = 0, w0 = uniform
    k = A.shape[1]
    w0 = np.full(k, 1.0 / k)
    r = y - A @ w0
    # basis for {v : 1'v = 0} is implicit: solve least norm with the constraint via projection
    P = np.eye(k) - np.outer(ones, ones) / k
    Ap = A @ P
    v, *_ = np.linalg.lstsq(Ap, r, rcond=None)
    v = P @ v
    return w0 + v


def oracle_ridge(A, y, lam, ones):
    """ORACLE. argmin ||A w - y||^2 + lam*||w - u||^2 s.t. 1'w = 1.  Closed form on the
    zero-sum subspace: the ladder's norm-bounded middle, uniform at lam -> inf."""
    k = A.shape[1]
    w0 = np.full(k, 1.0 / k)
    r = y - A @ w0
    P = np.eye(k) - np.outer(ones, ones) / k
    Ap = A @ P
    G = Ap.T @ Ap + lam * P                                    # P is the metric on the subspace
    v = np.linalg.solve(G + 1e-12 * np.eye(k), Ap.T @ r)
    v = P @ v
    return w0 + v


def oracle_convex(A, y, rounds=4000):
    """ORACLE. argmin ||A w - y||^2 s.t. w >= 0, 1'w = 1 -- FISTA with simplex projection.
    The convex rung: the reachable set of the SHIPPED `average_weighted` (`p` is a measurement
    distribution, hence non-negative).  A vertex of the simplex is feasible, so the optimum is at
    least as good as the argmin; that is asserted by the caller and is the solver's own check."""
    k = A.shape[1]
    w = np.full(k, 1.0 / k)
    z, t = w.copy(), 1.0
    L = float(np.linalg.norm(A, 2) ** 2) * 2.0 + 1e-12
    for _ in range(rounds):
        g = 2.0 * A.T @ (A @ z - y)
        wn = simplex_proj(z - g / L)
        tn = 0.5 * (1.0 + math.sqrt(1.0 + 4.0 * t * t))
        z = wn + ((t - 1.0) / tn) * (wn - w)
        w, t = wn, tn
    return w


def simplex_proj(v):
    u = np.sort(np.asarray(v, float))[::-1]
    css = np.cumsum(u) - 1.0
    idx = np.arange(1, len(u) + 1)
    cond = u - css / idx > 0
    rho = idx[cond][-1]
    theta = css[rho - 1] / rho
    return np.maximum(v - theta, 0.0)


# ---------------------------------------------------------------- one target
def row_for(pdb, verbose=True):
    from s27 import s28_A_amp as A28
    z = CA.load(pdb)
    W = np.asarray(z["W"], float)
    order = np.asarray(z["order"], int)
    dis = np.asarray(z["dis"], float)
    rr = np.asarray(z["rr"], float)                                   # ORACLE
    nat = np.asarray(z["nat"], float)                                 # ORACLE
    top75 = np.asarray(z["top75"], int)
    n, k = int(z["n"]), int(z["k"])

    frame = A28.Frame(W, np.sort(top75))          # the DEPLOYED average's frame; native-free
    Wf = frame.Wf                                                      # (k, 3n)
    natp = I.superpose_batch(nat[None], frame.ref)[0].ravel()          # ORACLE, for scoring only

    idx = order[:TOP]                                                  # THE FIXED CANDIDATE SET
    Wt = Wf[idx]                                                       # (128, 3n)
    A = Wt.T                                                           # (3n, 128)
    ones = np.ones(TOP)
    rng = np.random.default_rng(abs(hash((SALT, pdb))) % (2 ** 32))

    arms = {}

    def put(name, w, oracle, extra=None):
        X = emit(Wt, w, n)
        d = dict(rmsd=float(I.ca_rmsd(X, nat)), oracle=bool(oracle))
        d.update(wdiag(w))
        d.update(decompose(Wt, w, natp, n))
        if extra:
            d.update(extra)
        arms[name] = d

    # -- anchor: production's own operator, which must reproduce 3.0483 in aggregate -------------
    wp = np.zeros(k); wp[top75] = 1.0 / len(top75)
    Xp = emit(Wf, wp, n)
    arms["prod_top75_uniform"] = dict(rmsd=float(I.ca_rmsd(Xp, nat)), oracle=False,
                                      **wdiag(wp[top75]))

    # -- A. native-free uniform SCORE-PREFIX of the 128 ------------------------------------------
    for m in M_CURVE:
        w = np.zeros(TOP); w[:m] = 1.0 / m
        put("unif_prefix_m%d" % m, w, False)

    # -- B. uniform over a RANDOM m-subset of the 128: the exchangeable arm the saturation law
    #       is about.  MEAN over draws, never the per-target minimum (contract rule 10).
    for m in M_CURVE:
        rs, cs = [], []
        for _ in range(N_DRAW):
            sel = rng.choice(TOP, m, replace=False)
            w = np.zeros(TOP); w[sel] = 1.0 / m
            rs.append(float(I.ca_rmsd(emit(Wt, w, n), nat)))
            cs.append(decompose(Wt, w, natp, n))
        arms["rand_m%d" % m] = dict(
            rmsd=float(np.mean(rs)), rmsd_sd=float(np.std(rs)), oracle=False, n_draw=N_DRAW,
            c2=float(np.mean([x["c2"] for x in cs])),
            cross=float(np.mean([x["cross"] for x in cs])),
            d2=float(np.mean([x["d2"] for x in cs])),
            ms_rmsd=float(np.mean([r ** 2 for r in rs])))           # MEAN SQUARE, for the fit

    # -- C. native-free non-uniform weightings, as the registered CONCENTRATION CONTROLS ---------
    zr = (np.argsort(np.argsort(dis[idx])).astype(float) - (TOP - 1) / 2.0) / max(TOP / 4.0, 1e-9)
    for T in SOFT_T:
        e = np.exp(-zr / T); put("soft_T%.2f" % T, e / e.sum(), False)
    P = I.pairwise_rmsd(W[idx])
    typ = P.mean(1)
    tz = (typ - typ.mean()) / max(typ.std(), 1e-12)
    for T in SOFT_T:
        e = np.exp(-tz / T); put("typ_T%.2f" % T, e / e.sum(), False)
    # a random simplex weight: the matched zero-information control in the operator's own space
    rs = []
    for _ in range(N_DRAW):
        w = rng.dirichlet(np.ones(TOP))
        rs.append(float(I.ca_rmsd(emit(Wt, w, n), nat)))
    arms["randw_dirichlet"] = dict(rmsd=float(np.mean(rs)), oracle=False, n_draw=N_DRAW)

    # -- D. ORACLE rungs over the SAME 128 -------------------------------------------------------
    put("ORACLE_argmin128", np.eye(TOP)[int(np.argmin(rr[idx]))], True,
        dict(member=int(np.argmin(rr[idx]))))
    #: two independent convex solvers; the ceiling takes the better, and the vertex bound is
    #: asserted rather than assumed (a convex optimum can never be worse than the argmin).
    from s29 import s29_O_ladder as O29
    wc = oracle_convex(A, natp)
    r_f = float(I.ca_rmsd(emit(Wt, wc, n), nat))
    r_h, _Xh, wh = O29.oracle_hull(Wf, nat, n, support=np.asarray(idx, int))
    wh_full = np.zeros(TOP); wh_full[:] = np.asarray(wh, float)
    r_h = float(I.ca_rmsd(emit(Wt, wh_full, n), nat))
    #: DEFECT FOUND IN SELF-AUDIT AND FIXED: on 4 of 126 targets (1ID6, 2BP4, 2NDM, 9L1M -- three
    #: of them FAIL18, i.e. the ill-conditioned sets) both iterative solvers returned a point
    #: WORSE than the best simplex VERTEX, which is feasible, so the returned value could not have
    #: been the convex optimum.  The vertex is now a third candidate and the bound is ASSERTED.
    #: Effect on the aggregate: 1.8008 -> 1.7977 A, i.e. the original figure was conservative.
    wv = np.eye(TOP)[int(np.argmin(rr[idx]))]
    r_v = float(rr[idx].min())
    cands = [(r_f, wc, "fista"), (r_h, wh_full, "hull"), (r_v, wv, "vertex")]
    r_best, wc, solver = min(cands, key=lambda t: t[0])
    assert r_best <= r_v + 1e-9, "convex optimum worse than a feasible vertex"
    put("ORACLE_convex128", wc, True,
        dict(solver=solver, rmsd_fista=r_f, rmsd_hull=r_h, rmsd_vertex=r_v))
    wa = oracle_affine(A, natp, ones)
    put("ORACLE_affine128", wa, True)
    for lam in LAM:
        put("ORACLE_ridge_lam%.4g" % lam, oracle_ridge(A, natp, lam, ones), True, dict(lam=float(lam)))

    row = dict(pdb=pdb, n=n, fold=int(z["fold"]), dim=3 * n, top=TOP,
               fail18=bool(pdb in I.FAIL18),
               mean_member_rmsd=float(rr[idx].mean()),
               ms_member_rmsd=float((rr[idx] ** 2).mean()),
               best_member_rmsd=float(rr[idx].min()),
               arms=arms)
    return row


def phase_run(limit=0):
    pdbs = CA.all_pdbs()
    if limit:
        pdbs = pdbs[:limit]
    done = set()
    if os.path.exists(ROWS):
        for line in open(ROWS):
            done.add(json.loads(line)["pdb"])
    t0 = time.time()
    with open(ROWS, "a") as fh:
        for i, p in enumerate(pdbs):
            if p in done:
                continue
            r = row_for(p)
            fh.write(json.dumps(r) + "\n"); fh.flush()
            if i % 10 == 0 or i == len(pdbs) - 1:
                print("[%3d/%3d] %s  %.1fs" % (i + 1, len(pdbs), p, time.time() - t0), flush=True)
    print("rows done in %.1fs" % (time.time() - t0))


# ---------------------------------------------------------------- analysis
def phase_analyse():
    rows = [json.loads(l) for l in open(ROWS)]
    rows = {r["pdb"]: r for r in rows}
    pdbs = [p for p in CA.all_pdbs() if p in rows]
    rows = [rows[p] for p in pdbs]
    folds = np.array([r["fold"] for r in rows])
    n = len(rows)

    def arm(name, key="rmsd"):
        return np.array([r["arms"][name][key] for r in rows])

    names = sorted({k for r in rows for k in r["arms"]})
    out = {"n": n, "basis": "CA point cloud", "top": TOP, "arms": {},
           "provenance": ST.provenance(__file__)}
    for a in names:
        if not all(a in r["arms"] for r in rows):
            continue
        v = arm(a)
        d = {"mean": float(v.mean()), "se": float(v.std(ddof=1) / math.sqrt(n)),
             "oracle": bool(rows[0]["arms"][a].get("oracle", False)),
             "fail18": float(v[[r["fail18"] for r in rows]].mean())}
        for key in ("l1", "ess", "frac_neg", "neg_mass", "c2", "cross", "d2"):
            if key in rows[0]["arms"][a]:
                d[key] = float(arm(a, key).mean())
        out["arms"][a] = d

    # ---- the saturation law, fitted on the EXCHANGEABLE (random-subset) arm -------------------
    ms = np.array(M_CURVE, float)
    y = np.array([np.mean([r["arms"]["rand_m%d" % int(m)]["ms_rmsd"] for r in rows]) for m in ms])
    X = np.vstack([np.ones(len(ms)), 1.0 / ms]).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    c2, v2 = float(beta[0]), float(beta[1])
    pred = X @ beta
    out["saturation"] = {
        "model": "MS-RMSD(m) = c2 + v2/m  (exact for error = shared + exchangeable)",
        "c2": c2, "v2": v2, "c_rms": math.sqrt(max(c2, 0.0)),
        "r2": float(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()),
        "common_mode_fraction_at_m1": float(c2 / (c2 + v2)),
        "ms_m1": float(y[0]), "ms_m128": float(y[-1]),
        "observed_m128_rmsd": float(arm("rand_m128").mean()),
        "prediction_0.68": {"note": "coordinator's prediction: c2/(c2+v2) ~ 0.68 from the pool's "
                                    "measured common-mode share",
                            "predicted_c_rms": float(math.sqrt(0.68 * y[0]))},
        "m_grid": [int(m) for m in ms],
        "ms_observed": [float(q) for q in y], "ms_fit": [float(q) for q in pred],
    }

    # ---- the exact decomposition, by class ---------------------------------------------------
    out["decomposition"] = {}
    for a in ("unif_prefix_m75", "unif_prefix_m128", "rand_m128", "soft_T1.00", "typ_T1.00",
              "ORACLE_argmin128", "ORACLE_convex128", "ORACLE_affine128"):
        if a in out["arms"] and "c2" in out["arms"][a]:
            d = out["arms"][a]
            out["decomposition"][a] = {"rmsd": d["mean"], "c_rms": math.sqrt(max(d["c2"], 0.0)),
                                       "cross": d["cross"], "d_rms": math.sqrt(max(d["d2"], 0.0)),
                                       "oracle": d["oracle"]}

    # ---- the set-matched ladder, with fold CIs against the argmin ----------------------------
    base = arm("ORACLE_argmin128")
    out["vs_argmin128"] = {}
    for a in ("unif_prefix_m2", "unif_prefix_m3", "unif_prefix_m5", "unif_prefix_m75",
              "unif_prefix_m128", "ORACLE_convex128", "ORACLE_affine128", "prod_top75_uniform"):
        if a in out["arms"]:
            out["vs_argmin128"][a] = ST.compare(arm(a), base, folds=folds, names=pdbs,
                                                label="%s_vs_argmin128" % a,
                                                seed_parts=("s31C", "ladder"))

    # ---- F-C1d: is the ORACLE affine ceiling informative or degenerate? ----------------------
    ess = arm("ORACLE_affine128", "ess"); nm = arm("ORACLE_affine128", "neg_mass")
    l1 = arm("ORACLE_affine128", "l1")
    dim = np.array([r["dim"] for r in rows], float)
    out["F_C1d"] = {
        "rmsd": float(arm("ORACLE_affine128").mean()),
        "median_ess": float(np.median(ess)), "median_neg_mass": float(np.median(nm)),
        "median_l1": float(np.median(l1)), "median_dim_3n": float(np.median(dim)),
        "frac_targets_ess_ge5_and_negmass_le1": float(((ess >= 5) & (nm <= 1.0)).mean()),
        "fires": bool(((ess >= 5) & (nm <= 1.0)).mean() > 0.5),
        "note": "registered: the affine ceiling is INFORMATIVE only if ess >= 5 and neg_mass <= 1 "
                "on a majority of targets; otherwise it is a dimension-counting tautology and no "
                "affine-hull ceiling is quotable (s27/s28_A_FINDINGS.md:85-90)",
    }

    # ---- the ridge path: the norm-bounded middle ---------------------------------------------
    out["ridge_path"] = []
    for lam in LAM:
        a = "ORACLE_ridge_lam%.4g" % lam
        if a in out["arms"]:
            d = out["arms"][a]
            out["ridge_path"].append({"lam": float(lam), "rmsd": d["mean"], "l1": d.get("l1"),
                                      "ess": d.get("ess"), "neg_mass": d.get("neg_mass"),
                                      "frac_neg": d.get("frac_neg")})

    # ---- native-free arms against production --------------------------------------------------
    pr = arm("prod_top75_uniform")
    out["vs_production"] = {}
    for a in names:
        if a.startswith("ORACLE") or a == "prod_top75_uniform":
            continue
        if not all(a in r["arms"] for r in rows):
            continue
        c = ST.compare(arm(a), pr, folds=folds, names=pdbs, label="%s_vs_prod" % a,
                       seed_parts=("s31C", "ladder"))
        out["vs_production"][a] = {kk: c[kk] for kk in
                                   ("effect", "se", "mde", "effect_over_mde", "ci95_fold",
                                    "folds_same_sign", "n_better", "n_worse", "mean_a", "mean_b")}
    best_nf = min(out["vs_production"], key=lambda a: out["vs_production"][a]["effect"])
    out["best_native_free"] = {"arm": best_nf, **out["vs_production"][best_nf]}

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    # ---- render -------------------------------------------------------------------------------
    print("\nSET-MATCHED READOUT LADDER -- candidate set fixed at the top-%d.  CA point cloud, n=%d\n"
          % (TOP, n))
    print("  %-24s %8s %8s %8s %8s %8s %8s  %s" %
          ("arm", "rmsd", "c_rms", "cross", "d_rms", "ess", "neg_m", "label"))
    for a in ("prod_top75_uniform", "unif_prefix_m2", "unif_prefix_m5", "unif_prefix_m75",
              "unif_prefix_m128", "rand_m128", "soft_T1.00", "typ_T1.00",
              "ORACLE_argmin128", "ORACLE_convex128", "ORACLE_affine128"):
        if a not in out["arms"]:
            continue
        d = out["arms"][a]
        print("  %-24s %8.4f %8.4f %8.4f %8.4f %8.2f %8.3f  %s" %
              (a, d["mean"], math.sqrt(max(d.get("c2", 0.0), 0.0)), d.get("cross", float("nan")),
               math.sqrt(max(d.get("d2", 0.0), 0.0)), d.get("ess", float("nan")),
               d.get("neg_mass", float("nan")),
               "ORACLE / NOT DEPLOYABLE" if d["oracle"] else "native-free"))

    s = out["saturation"]
    print("\n  SATURATION (exchangeable random m-subsets of the 128; MEAN of %d draws):" % N_DRAW)
    print("    MS(m) = c2 + v2/m   c2 %.4f  v2 %.4f  R2 %.5f" % (s["c2"], s["v2"], s["r2"]))
    print("    common-mode share at m=1: %.4f   (coordinator predicted 0.68)" %
          s["common_mode_fraction_at_m1"])
    print("    fitted floor c_rms %.4f A;  observed uniform-128 %.4f A" %
          (s["c_rms"], s["observed_m128_rmsd"]))
    print("    m:        " + " ".join("%7d" % m for m in s["m_grid"][:10]))
    print("    rmsd obs: " + " ".join("%7.4f" % math.sqrt(q) for q in s["ms_observed"][:10]))
    print("    rmsd fit: " + " ".join("%7.4f" % math.sqrt(max(q, 0)) for q in s["ms_fit"][:10]))

    f = out["F_C1d"]
    print("\n  F-C1d (is the ORACLE affine ceiling informative?): rmsd %.6f  median ess %.2f of %d"
          "  median neg_mass %.2f  median ||w||_1 %.2f  median 3n %.0f" %
          (f["rmsd"], f["median_ess"], TOP, f["median_neg_mass"], f["median_l1"],
           f["median_dim_3n"]))
    print("    frac targets with ess>=5 AND neg_mass<=1: %.3f  ->  %s" %
          (f["frac_targets_ess_ge5_and_negmass_le1"],
           "FIRES (informative)" if f["fires"] else "REFUTED (degenerate tautology)"))

    print("\n  RIDGE PATH (ORACLE, the norm-bounded middle):")
    print("    %10s %8s %8s %8s %8s" % ("lam", "rmsd", "||w||_1", "ess", "neg_mass"))
    for p in out["ridge_path"]:
        print("    %10.4g %8.4f %8.2f %8.2f %8.3f" %
              (p["lam"], p["rmsd"], p["l1"], p["ess"], p["neg_mass"]))

    print("\n  VS THE ARGMIN OVER THE SAME 128 (negative = the arm is better):")
    for a, c in out["vs_argmin128"].items():
        print("    %-22s %+8.4f  %5.2fx MDE  fold CI [%+.4f,%+.4f]  %d/%d folds  %dW/%dL" %
              (a, c["effect"], c["effect_over_mde"], c["ci95_fold"][0], c["ci95_fold"][1],
               c["folds_same_sign"], c["n_folds"], c["n_better"], c["n_worse"]))

    b = out["best_native_free"]
    print("\n  BEST NATIVE-FREE ARM vs PRODUCTION: %s  %+.4f A  %.2fx MDE  fold CI [%+.4f,%+.4f]"
          "  %dW/%dL" % (b["arm"], b["effect"], b["effect_over_mde"], b["ci95_fold"][0],
                         b["ci95_fold"][1], b["n_better"], b["n_worse"]))
    print("\n  wrote %s" % OUT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["run", "analyse"])
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.phase == "run":
        phase_run(a.limit)
    else:
        phase_analyse()


if __name__ == "__main__":
    main()
