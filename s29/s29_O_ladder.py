#!/usr/bin/env python
"""s29/s29_O_ladder.py -- S29 lane O: the ORACLE ceiling ladder by operator class, and the
typicality-axis probe (rung 6).  Pre-registration: `s29/PREREG_S29_O.md` (read it first).

EVERY NUMBER HERE IS ORACLE (it reads `nat_ca` / `oracle_rr`) EXCEPT the leave-fold-out step
of rung 6 (`lfo_*`), whose parameter is chosen on the other four folds' natives and applied
to the fifth: that arm is DEPLOYABLE and is the lane's one paired contrast.  ORACLE reads
happen only inside functions whose names begin `oracle_` or in `cloud_row` under the
`# ORACLE` markers; the NaN-poison test in `tests/test_s29_O.py` asserts the LFO structure
is bit-identical when the held-out target's native is NaN.

The ladder (point cloud AND built chain, 126 targets):
  1. best single member (top-75 / pool)                        items best1_top75, best1_pool
  2. best m-subset of the DIS order, m = 1..500                item  bestm
  3. best basin average, k in {2,3,4,6,8}, top-75 / top-500     items basin_<S>_k<k>   (10)
  4. best sparse convex combination, s in {2,3,5,10,20}         items sparse_<S>_s<s>  (10)
  5. best convex combination (the hull), top-75 / 500           items hull_top75, hull_pool
  6. the typicality axis: u = prod - blind, v = native - prod;  items prod, lfo_LIB75, lfo_BPRIME
     cos(u, v), random-field reference, ORACLE global / per-target step, LEAVE-FOLD-OUT step

Phases (each a governed job, per-target / per-item checkpoints, resumable):
    python s29/s29_O_ladder.py cloud   [--pdbs A,B] [--limit N]   # point cloud, structures stored
    python s29/s29_O_ladder.py lfo                                # rung 6 fold choices (needs 126 rows)
    python s29/s29_O_ladder.py chain   [--groups A,B,C,D]         # built chain of stored structures
    python s29/s29_O_ladder.py analyse                            # summary json + ST.fmt blocks + table
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

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

RESULTS = os.path.join(HERE, "results")
STRUCTS = os.path.join(RESULTS, "s29_O_structs")
CLOUD_ROWS = os.path.join(RESULTS, "s29_O_cloud_rows.jsonl")
CHAIN_ROWS = os.path.join(RESULTS, "s29_O_chain_rows.jsonl")
LFO_JSON = os.path.join(RESULTS, "s29_O_lfo.json")
SUMMARY_JSON = os.path.join(RESULTS, "s29_O_summary.json")
SALT = "s29O"
K, M = 500, 75
BASIN_K = (2, 3, 4, 6, 8)
SPARSE_S = (2, 3, 5, 10, 20)
SETS = ("top75", "pool")
T_GRID = np.round(np.arange(-1.0, 2.0 + 1e-9, 0.1), 6)        # 31 values, fixed in the prereg
N_RAND_FIELDS = 16
BLIND_DEFS = ("LIB75", "BPRIME")
#: rung 8 (lane T's S29-L11 prediction 4): the one-parameter family X(eta) = production + eta*PC1,
#: PC1 the top shape mode of the top-75 members' deviations from the production average, scaled so
#: that eta is in ANGSTROM of point-cloud RMSD displacement.  Grid fixed here, before any number.
ETA_GRID = np.round(np.arange(-3.0, 3.0 + 1e-9, 0.1), 6)
PC1_ROWS = os.path.join(RESULTS, "s29_O_pc1_rows.jsonl")
PC1_JSON = os.path.join(RESULTS, "s29_O_pc1.json")
#: rung 9 (lane M's C13): the deployed quantum stage widens the retained prefix to 2**n = 128
#: (`core/pipeline.py:758`), and the set-equality theorem (S25) says the CVaR tail's support is
#: always a PREFIX of the energy order.  So the ORACLE ceiling of the top-128 prefix bounds every
#: quantum arm this project has run or could run in the deployed encoding.
M128 = 128
P128_ROWS = os.path.join(RESULTS, "s29_O_p128_rows.jsonl")
P128_JSON = os.path.join(RESULTS, "s29_O_p128.json")
RHO_SUM = 1e3                                                  # the sum-to-one augmented row
HULL_ROUNDS, SPARSE_ROUNDS = 5, 3
CHAIN_GROUPS = {
    "A": ["prod", "lfo_LIB75", "lfo_BPRIME"],
    "E": ["best1_top128", "bestm128", "hull_top128"],
    "B": ["best1_top75", "best1_pool", "bestm", "hull_top75", "hull_pool"],
    "C": ["sparse_%s_s%d" % (S, s) for S in SETS for s in SPARSE_S],
    "D": ["basin_%s_k%d" % (S, k) for S in SETS for k in BASIN_K],
}
ALL_ITEMS = [it for g in "ABCDE" for it in CHAIN_GROUPS[g]]


# ================================================================== pool, frames, geometry
def load_pool(pdb):
    """The shipped pool, DIS, the production top-75 (tie key) and the full DIS order."""
    from s27 import run_pool as RP
    from s27 import s28_A_amp as A
    cand, dis, top, dg = A.load_pool(pdb)
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    order = np.lexsort((key, dis)).astype(int)
    assert set(order[:M].tolist()) == set(top.tolist()), "DIS order head != production top-75"
    return cand, dis, top, order, dg


def superpose_one(X, ref):
    return I.superpose_batch(np.asarray(X, float)[None], np.asarray(ref, float))[0]


def avg_of(W, idx, P):
    """`s12.instrument.coordinate_average` on W[idx] using the precomputed pairwise matrix P."""
    idx = np.asarray(idx, int)
    sub = P[np.ix_(idx, idx)]
    b = idx[int(np.argmin(sub.mean(1)))]
    return I.superpose_batch(W[idx], W[b]).mean(0), int(b)


def rigid_basis(C):
    """Orthonormal (3n, 6) basis of the rigid-body tangent space at C (3 translations, 3
    infinitesimal rotations about C's centroid)."""
    C = np.asarray(C, float)
    n = len(C)
    Cc = C - C.mean(0)
    cols = []
    for a in range(3):
        e = np.zeros((n, 3)); e[:, a] = 1.0
        cols.append(e.ravel())
    for a in range(3):
        w = np.zeros(3); w[a] = 1.0
        cols.append(np.cross(w[None, :], Cc).ravel())
    B = np.column_stack(cols)
    Q, _ = np.linalg.qr(B)
    return Q


def remove_rigid(d, C):
    """Project a displacement field d (3n,) onto the complement of the rigid-body tangent space."""
    Q = rigid_basis(C)
    d = np.asarray(d, float).ravel()
    return d - Q @ (Q.T @ d)


def cosine(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else float("nan")


def struct_diag(C):
    C = np.asarray(C, float)
    if not np.isfinite(C).all():
        return dict(rg=float("nan"), bond=float("nan"))
    return dict(rg=float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean())),
                bond=float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean()))


# ============================================================ convex solvers (rungs 4, 5)
def convex_nnls(A, y, rho=RHO_SUM):
    """min |A w - y| over w >= 0, sum w = 1 (the sum as a weighted augmented row; w renormalised
    exactly).  A is (3n, s) columns = flattened posed windows, y the posed target (3n,)."""
    from scipy.optimize import nnls
    A = np.asarray(A, float); y = np.asarray(y, float).ravel()
    Aa = np.vstack([A, rho * np.ones((1, A.shape[1]))])
    ya = np.concatenate([y, [rho]])
    w, _ = nnls(Aa, ya, maxiter=50 * max(A.shape))
    s = w.sum()
    if s <= 0:                                             # degenerate: fall back to uniform
        w = np.full(A.shape[1], 1.0 / A.shape[1])
    else:
        w = w / s
    return w


def oracle_hull(Wf, nat, n, rounds=HULL_ROUNDS, support=None, natp0=None):
    """ORACLE.  Best convex combination of the posed windows Wf (k, 3n) against the native,
    the native's pose re-solved on the current combination each round (S10-5's joint `cf`
    problem by alternation).  Returns (rmsd after Kabsch, X (n,3), w over `support`)."""
    Wf = np.asarray(Wf, float)
    idx = np.arange(len(Wf)) if support is None else np.asarray(support, int)
    A = Wf[idx].T                                            # (3n, s)
    natp = (superpose_one(nat, Wf[idx[0]].reshape(n, 3)) if natp0 is None else natp0).ravel()
    w = None
    for _ in range(int(rounds)):
        w = convex_nnls(A, natp)
        X = (w @ Wf[idx]).reshape(n, 3)
        natp = superpose_one(nat, X).ravel()
    X = (w @ Wf[idx]).reshape(n, 3)
    return float(I.ca_rmsd(X, nat)), X, w


def oracle_sparse_greedy(Wf, nat, n, s_list=SPARSE_S, rounds=SPARSE_ROUNDS, natp0=None):
    """ORACLE.  Greedy forward support selection with the convex NNLS refit; nested path, so one
    pass returns the structure at every s in `s_list`.  Returns {s: (rmsd, X, support, w)}."""
    Wf = np.asarray(Wf, float)
    k = len(Wf)
    s_max = int(max(s_list))
    natp = (superpose_one(nat, Wf[0].reshape(n, 3)) if natp0 is None else natp0).ravel()
    support: List[int] = []
    out = {}
    w = None
    for step in range(1, min(s_max, k) + 1):
        best = None
        for j in range(k):
            if j in support:
                continue
            cols = support + [j]
            wj = convex_nnls(Wf[cols].T, natp)
            res = float(np.linalg.norm(wj @ Wf[cols] - natp))
            if best is None or res < best[0] - 1e-12:
                best = (res, j, wj)
        support.append(best[1]); w = best[2]
        for _ in range(int(rounds)):                         # frame alternation on the new support
            X = (w @ Wf[support]).reshape(n, 3)
            natp = superpose_one(nat, X).ravel()
            w = convex_nnls(Wf[support].T, natp)
        if step in s_list:
            X = (w @ Wf[support]).reshape(n, 3)
            out[step] = (float(I.ca_rmsd(X, nat)), X, list(support), w.copy())
    return out


# ======================================================================= rungs 2 and 3
def oracle_m_curve(W, order, P, nat, m_max=K):
    """ORACLE.  r(m) for the uniform medoid-frame average of the first m of `order`, m = 1..m_max,
    and the structures at every m (returned lazily: only the curve; call avg_of for a chosen m)."""
    r = np.empty(int(m_max))
    for m in range(1, int(m_max) + 1):
        X, _ = avg_of(W, order[:m], P)
        r[m - 1] = I.ca_rmsd(X, nat)
    return r


def cluster_labels(P_sub, k):
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    D = np.asarray(P_sub, float)
    D = 0.5 * (D + D.T); np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    return fcluster(Z, int(k), criterion="maxclust")


def oracle_best_basin(W, idx, P, nat, k):
    """ORACLE.  Cluster W[idx] into k basins (average linkage on pairwise CA-RMSD); average every
    basin of size >= 2; return the best (rmsd, X, size, n_basins_found, n_eligible)."""
    idx = np.asarray(idx, int)
    lab = cluster_labels(P[np.ix_(idx, idx)], k)
    best = None
    n_found = len(set(lab.tolist())); n_elig = 0
    for l in sorted(set(lab.tolist())):
        mem = idx[lab == l]
        if len(mem) < 2:
            continue
        n_elig += 1
        X, _ = avg_of(W, mem, P)
        r = float(I.ca_rmsd(X, nat))
        if best is None or r < best[0]:
            best = (r, X, int(len(mem)))
    if best is None:
        return float("nan"), np.full((W.shape[1], 3), np.nan), 0, n_found, 0
    return best[0], best[1], best[2], n_found, n_elig


# ================================================================ rung 6: the blind clouds
def _bias_native_frame(C, nat):
    """S24's `_bias`: the average posed on the native, minus the native (native frame)."""
    C = np.asarray(C, float); nat = np.asarray(nat, float)
    Cc = C - C.mean(0); Nc = nat - nat.mean(0)
    U, S, Vt = np.linalg.svd(Cc.T @ Nc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return Cc @ R.T - Nc


def blind_cloud(pdb, u, dg, which):
    """Rebuild S24's blind average exactly as S24 built it (its own stable RNG).  Returns
    (B (n,3), members (75,) universe indices, meta)."""
    Wall = np.asarray(u["W"], float); N = len(Wall); n = int(u["n"])
    if which == "LIB75":
        rng = SD.stable_rng("biasalign", pdb)
        c1 = rng.choice(N, M, replace=False)                 # biasalign.py's C1, the first draw
        B, _ = I.coordinate_average(Wall[c1])
        return B, c1, dict(source="s24/biasalign.py C1 (uniform 75 of the universe)")
    if which == "BPRIME":
        idx = I.pool_idx(u)
        i, j = I.pair_index(n)
        rng = SD.stable_rng("qmatch", pdb)
        mask = np.ones(N, bool); mask[idx] = False
        pool_free = np.flatnonzero(mask)
        pick = rng.choice(pool_free, min(2000, len(pool_free)), replace=False)
        scf = np.asarray(I.shipped_score(dg, I.pair_dists(Wall[pick], i, j)), float)
        mem = pick[np.argsort(scf, kind="stable")[:M]]
        B, _ = I.coordinate_average(Wall[mem])
        return B, mem, dict(source="s24/qmatch.py B' (top-75 by the shipped score of 2000 uniform, pool excluded)")
    raise ValueError(which)


def s24_gate_values(pdb, u, dg, B, which):
    """The S24 artefact values this cloud must reproduce, and the values it produces."""
    Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float); n = int(u["n"])
    i, j = I.pair_index(n)
    idx = I.pool_idx(u); Wpool = Wall[idx]
    sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
    A24 = Wpool[np.argsort(sc, kind="stable")[:M]]                # S24's A (stable argsort)
    CA24, _ = I.coordinate_average(A24)
    eA = _bias_native_frame(CA24, nat); eB = _bias_native_frame(B, nat)
    got = dict(rmsd_B=float(I.ca_rmsd(B, nat)), cos=cosine(eA, eB), rmsd_A=float(I.ca_rmsd(CA24, nat)))
    if which == "LIB75":
        f = os.path.join(ROOT, "s24", "results", "biasalign.json")
        row = {r["pdb"]: r for r in json.load(open(f))["rows"]}[pdb]
        want = dict(rmsd_B=float(row["mix"]["m0_75"]), cos=float(row["cos_A_C"]), rmsd_A=float(row["rmsd_A"]))
    else:
        f = os.path.join(ROOT, "s24", "results", "qmatch.json")
        row = {r["pdb"]: r for r in json.load(open(f))["rows"]}[pdb]
        want = dict(rmsd_B=float(row["rmsd_B"]), cos=float(row["cos_A_B"]), rmsd_A=float(row["rmsd_A"]))
    return got, want, CA24


def axis_probe(C, B, nat, pdb, which, t_grid=T_GRID, n_fields=N_RAND_FIELDS):
    """Rung 6 per target.  C production cloud (its own frame); B the blind average; nat ORACLE.
    Returns the per-target record and the rigid-removed u (3n,) for the step structures."""
    C = np.asarray(C, float); n = len(C)
    B_sup = superpose_one(B, C)
    N_sup = superpose_one(nat, C)                                             # ORACLE
    u = remove_rigid((C - B_sup).ravel(), C)
    v = remove_rigid((N_sup - C).ravel(), C)                                  # ORACLE
    cos_uv = cosine(u, v)                                                     # ORACLE
    rng = SD.stable_rng(pdb, "s29O_randfield", which, salt=SALT)
    cos_rand = []
    for _ in range(int(n_fields)):
        g = remove_rigid(rng.normal(0.0, 1.0, 3 * n), C)
        cos_rand.append(cosine(g, v))                                         # ORACLE
    cos_rand = np.asarray(cos_rand)
    curve = np.array([I.ca_rmsd(C + t * u.reshape(n, 3), nat) for t in t_grid])   # ORACLE
    rec = dict(cos_uv=cos_uv, cos_rand_mean=float(cos_rand.mean()), cos_rand_absmean=float(np.abs(cos_rand).mean()),
               cos_rand_analytic_absmean=float(math.sqrt(2.0 / (math.pi * (3 * n - 6)))),
               norm_u=float(np.linalg.norm(u)), norm_v=float(np.linalg.norm(v)),
               rmsd_blind=float(I.ca_rmsd(B, nat)), rmsd_blind_in_C_frame_at_t_minus1=float(curve[0]),
               curve=curve.tolist(), t_grid=[float(t) for t in t_grid],
               t_best=float(t_grid[int(np.argmin(curve))]), r_best=float(curve.min()),
               r_t0=float(curve[int(np.argmin(np.abs(t_grid)))]))
    return rec, u


def lfo_structure(C, u, t):
    """The DEPLOYABLE step structure: production cloud plus t times the native-free axis u.
    Reads no native.  `t` comes from other folds."""
    C = np.asarray(C, float)
    return C + float(t) * np.asarray(u, float).reshape(len(C), 3)


# ============================================================== rung 8: the PC1 family
def pool_pc1(W, top, P, C):
    """NATIVE-FREE.  The pool's first shape mode: the top right-singular vector of the top-75
    members' deviations from the production average, in the average's own (medoid) frame, with the
    rigid-body component removed and scaled so that ||PC1||_2 = sqrt(n) -- i.e. adding eta*PC1
    moves the point cloud by exactly |eta| A of RMSD.  The SIGN is fixed by a deterministic,
    native-free convention (the largest-|component| entry is positive): lane T's point is that the
    marginals cannot supply the sign, so an achievable arm must carry an arbitrary one."""
    top = np.asarray(top, int)
    n = W.shape[1]
    idx = top[int(np.argmin(P[np.ix_(top, top)].mean(1)))]
    Wp = I.superpose_batch(W[top], W[idx])                     # the deployed frame
    D = (Wp - C[None]).reshape(len(top), 3 * n)
    U, S, Vt = np.linalg.svd(D - D.mean(0, keepdims=True), full_matrices=False)
    v = remove_rigid(Vt[0], C)
    nv = np.linalg.norm(v)
    if nv <= 0:
        return np.zeros(3 * n), dict(sd_along=0.0, var_frac=0.0, sv=[])
    v = v / nv
    v = v * (1.0 if v[int(np.argmax(np.abs(v)))] >= 0 else -1.0)
    proj = (D - D.mean(0, keepdims=True)) @ v
    var_frac = float(S[0] ** 2 / max((S ** 2).sum(), 1e-300))
    v = v * math.sqrt(n)                                       # eta in A of RMSD
    return v, dict(sd_along=float(proj.std(ddof=1) / math.sqrt(n)), var_frac=var_frac,
                   sv=[float(x) for x in S[:5]])


def pc1_row(pdb, verbose=True):
    """Rung 8 per target.  PC1 is NATIVE-FREE; the curve over eta is ORACLE."""
    t0 = time.time()
    cand, dis, top, order, dg = load_pool(pdb)
    W = cand.W; n = cand.n; nat = cand.nat_ca                                 # ORACLE label
    P = I.pairwise_rmsd(W)
    C, _ = avg_of(W, np.sort(top), P)
    v, meta = pool_pc1(W, np.sort(top), P, C)
    curve = np.array([I.ca_rmsd(C + e * v.reshape(n, 3), nat) for e in ETA_GRID])   # ORACLE
    j0 = int(np.argmin(np.abs(ETA_GRID)))
    row = dict(pdb=pdb, n=n, fold=int(cand.fold), fail18=bool(pdb in I.FAIL18),
               curve=curve.tolist(), eta_grid=[float(e) for e in ETA_GRID],
               r_prod=float(curve[j0]), eta_best=float(ETA_GRID[int(np.argmin(curve))]),
               r_best=float(curve.min()), secs=time.time() - t0, **meta)
    if verbose:
        print("  %s n=%d prod %.3f  PC1 sd %.3f (var frac %.3f)  ORACLE best eta %+.1f -> %.3f  %.1fs"
              % (pdb, n, row["r_prod"], meta["sd_along"], meta["var_frac"], row["eta_best"],
                 row["r_best"], row["secs"]), flush=True)
    return row


def phase_pc1(pdbs):
    done = jsonl_rows(PC1_ROWS)
    for pdb in pdbs:
        if (pdb,) in done:
            continue
        append_row(PC1_ROWS, pc1_row(pdb))
    rows = jsonl_rows(PC1_ROWS)
    if any((p,) not in rows for p in all_pdbs()):
        print("pc1 rows:", PC1_ROWS, "(partial)")
        return None
    return analyse_pc1()


def analyse_pc1():
    rows = jsonl_rows(PC1_ROWS)
    pdbs = all_pdbs()
    folds = ST.pinned_folds(pdbs)
    R = np.array([rows[(p,)]["curve"] for p in pdbs])                          # ORACLE curves
    prod = np.array([rows[(p,)]["r_prod"] for p in pdbs])
    j0 = int(np.argmin(np.abs(ETA_GRID)))
    assert np.allclose(R[:, j0], prod), "eta = 0 column != production"
    mean_curve = R.mean(0)
    jg = int(np.argmin(mean_curve))
    t_of, t_fold, ties = lfo_choices(R, folds, ETA_GRID)
    arm = np.array([R[i, int(np.argmin(np.abs(ETA_GRID - t_of[i])))] for i in range(len(pdbs))])
    pos = R[:, ETA_GRID >= 0].min(1); neg = R[:, ETA_GRID <= 0].min(1)
    eta_best = np.array([rows[(p,)]["eta_best"] for p in pdbs])
    out = dict(
        eta_grid=[float(e) for e in ETA_GRID], mean_curve=mean_curve.tolist(),
        oracle_global_eta=float(ETA_GRID[jg]), oracle_global_mean=float(mean_curve[jg]),
        oracle_global_gain=float(mean_curve[jg] - prod.mean()),
        oracle_per_target_mean=float(R.min(1).mean()),
        oracle_per_target_gain=float(R.min(1).mean() - prod.mean()),
        oracle_per_target_bok=ST.best_of_k_within(R, seed_parts=("s29O", "pc1eta")),
        oracle_sign_plus_mean=float(pos.mean()), oracle_sign_minus_mean=float(neg.mean()),
        frac_eta_best_positive=float((eta_best > 0).mean()),
        median_abs_eta_best=float(np.median(np.abs(eta_best))),
        eta_fold=t_fold, ties=ties, lfo_mean=float(arm.mean()), prod_mean=float(prod.mean()),
        lfo_cloud=strata(arm, prod, pdbs, folds,
                         "PC1 one-parameter family, LEAVE-FOLD-OUT eta vs prod (POINT CLOUD, DEPLOYABLE)"),
        oracle_global_cloud=strata(R[:, jg], prod, pdbs, folds,
                                   "ORACLE PC1 global eta vs prod (POINT CLOUD)"),
        oracle_per_target_cloud=strata(R.min(1), prod, pdbs, folds,
                                       "ORACLE PC1 per-target eta vs prod (POINT CLOUD, an order statistic)"),
        sd_along_mean=float(np.mean([rows[(p,)]["sd_along"] for p in pdbs])),
        var_frac_mean=float(np.mean([rows[(p,)]["var_frac"] for p in pdbs])),
        arm_cloud={p: float(v) for p, v in zip(pdbs, arm)},
        provenance=ST.provenance(__file__))
    ST.save_atomic(PC1_JSON, out, module_file=__file__)
    print("== rung 8 (PC1, lane T's one-parameter family). PC1 is NATIVE-FREE; every eta below is "
          "ORACLE except the leave-fold-out arm.")
    print("   pool sd along PC1 %.3f A, PC1 variance share %.3f"
          % (out["sd_along_mean"], out["var_frac_mean"]))
    print("   ORACLE global eta %+.1f -> %.4f (%+.4f vs production %.4f); ORACLE per-target eta %.4f "
          "(%+.4f, an order statistic); ORACLE best eta positive on %.0f%% of targets, median |eta| %.1f"
          % (out["oracle_global_eta"], out["oracle_global_mean"], out["oracle_global_gain"],
             out["prod_mean"], out["oracle_per_target_mean"], out["oracle_per_target_gain"],
             100 * out["frac_eta_best_positive"], out["median_abs_eta_best"]))
    print("   ORACLE best over eta >= 0 only %.4f; over eta <= 0 only %.4f (the one-global-sign arms)"
          % (out["oracle_sign_plus_mean"], out["oracle_sign_minus_mean"]))
    print("   per-target eta pricing: observed %+.4f, across-target null %+.4f (%.0f%% accounted), "
          "split-half %+.4f (%.0f%%), k_eff %.1f"
          % (out["oracle_per_target_bok"]["observed_gain"], out["oracle_per_target_bok"]["null_across_targets"],
             100 * out["oracle_per_target_bok"]["share_accounted"], out["oracle_per_target_bok"]["split_half"],
             100 * out["oracle_per_target_bok"]["split_half_frac"], out["oracle_per_target_bok"]["k_eff"]))
    print("   leave-fold-out eta per fold %s (ties %s)" % (out["eta_fold"], out["ties"]))
    for k in ("oracle_global_cloud", "oracle_per_target_cloud", "lfo_cloud"):
        print(ST.fmt(out[k]["all"]))
    print("pc1:", PC1_JSON)
    return out


# =========================================== rung 9: the top-128 prefix (the quantum field of view)
def p128_row(pdb, verbose=True):
    """Rung 9 per target. The prefix is NATIVE-FREE (the DIS order with the stable tie key);
    every RMSD below is ORACLE."""
    t0 = time.time()
    cand, dis, top, order, dg = load_pool(pdb)
    W = cand.W; n = cand.n; nat = cand.nat_ca; rr = cand.oracle_rr            # ORACLE labels
    P = I.pairwise_rmsd(W)
    pre = order[:M128]
    structs = {}
    # (a) best single member of the prefix                                    # ORACLE
    r = rr[pre]; m = r.min(); tie = np.isclose(r, m, atol=1e-9)
    j = int(pre[np.flatnonzero(tie)[0]])
    structs["best1_top128"] = W[j]
    # (b) best PREFIX-m average inside the top-128: the exact reachable-set ceiling of the CVaR
    #     tail, whose support is always a prefix of the energy order (S25's set-equality theorem)
    curve = np.array([I.ca_rmsd(avg_of(W, order[:mm], P)[0], nat) for mm in range(1, M128 + 1)])
    m_best = int(np.argmin(curve)) + 1
    structs["bestm128"] = avg_of(W, order[:m_best], P)[0]
    # (c) the convex hull of the prefix: upper-bounds ANY weighting inside the field of view
    from s27 import s28_A_amp as A
    frame = A.Frame(W, np.sort(top))
    natp0 = superpose_one(nat, frame.ref)
    r_h, X_h, w_h = oracle_hull(frame.Wf[np.sort(pre)], nat, n, natp0=natp0)   # ORACLE
    structs["hull_top128"] = X_h
    save_structs(pdb, structs, cand.seq, cand.fold, cand.n)
    row = dict(pdb=pdb, n=n, fold=int(cand.fold), fail18=bool(pdb in I.FAIL18),
               best1_top128=float(I.ca_rmsd(W[j], nat)), best1_member=j, best1_n_tied=int(tie.sum()),
               bestm128=float(curve[m_best - 1]), m128=m_best, curve128=curve.tolist(),
               hull_top128=r_h, hull_n_support=int((w_h > 1e-9).sum()),
               prod=float(curve[M - 1] if M <= M128 else np.nan),
               overlap_top75_in_128=float(len(set(pre.tolist()) & set(top.tolist())) / float(M)),
               secs=time.time() - t0)
    if verbose:
        print("  %s n=%d prod %.3f | top-128 best1 %.3f  best prefix-m %.3f (m=%d)  hull %.3f (supp %d)  %.1fs"
              % (pdb, n, row["prod"], row["best1_top128"], row["bestm128"], row["m128"],
                 row["hull_top128"], row["hull_n_support"], row["secs"]), flush=True)
    return row


def phase_p128(pdbs):
    done = jsonl_rows(P128_ROWS)
    for pdb in pdbs:
        if (pdb,) in done:
            continue
        append_row(P128_ROWS, p128_row(pdb))
    rows = jsonl_rows(P128_ROWS)
    if any((p,) not in rows for p in all_pdbs()):
        print("p128 rows:", P128_ROWS, "(partial)")
        return None
    return analyse_p128()


def analyse_p128():
    rows = jsonl_rows(P128_ROWS)
    cloud = jsonl_rows(CLOUD_ROWS)
    pdbs = all_pdbs()
    folds = ST.pinned_folds(pdbs)
    fm = np.array([p in I.FAIL18 for p in pdbs])
    G = lambda k: np.array([rows[(p,)][k] for p in pdbs])
    C = lambda it: np.array([cloud[(p,)]["items"][it]["rmsd_cloud"] for p in pdbs])
    prod = G("prod")
    Mc = np.array([cloud[(p,)]["m_curve"] for p in pdbs])                      # ORACLE, m = 1..500
    out = dict(prod_mean=float(prod.mean()), provenance=ST.provenance(__file__),
               overlap_mean=float(G("overlap_top75_in_128").mean()))
    trip = {}
    for nm, v in (("best1_top75", C("best1_top75")), ("best1_top128", G("best1_top128")),
                  ("best1_pool", C("best1_pool")),
                  ("bestm_top75", Mc[:, :M].min(1)), ("bestm128", G("bestm128")),
                  ("bestm_pool", Mc.min(1)),
                  ("hull_top75", C("hull_top75")), ("hull_top128", G("hull_top128")),
                  ("hull_pool", C("hull_pool"))):
        trip[nm] = dict(mean=float(v.mean()), median=float(np.median(v)), fail18=float(v[fm].mean()),
                        other108=float(v[~fm].mean()), frac_under2=float((v < 2.0).mean()),
                        vs_prod=ST.compare(v, prod, folds, names=pdbs,
                                           label="ORACLE %s vs production (POINT CLOUD)" % nm))
    out["arms"] = trip
    out["m128_median"] = float(np.median(G("m128")))
    out["m128_quartiles"] = [float(q) for q in np.percentile(G("m128"), [10, 25, 50, 75, 90])]
    out["frac_m128_at_128"] = float((G("m128") == M128).mean())
    out["hull_support_median"] = float(np.median(G("hull_n_support")))
    #: the three prefixes read against each other (the coordinator's question)
    for a, b in (("best1_top128", "best1_top75"), ("bestm128", "bestm_top75"),
                 ("hull_top128", "hull_top75"), ("best1_pool", "best1_top128"),
                 ("bestm_pool", "bestm128"), ("hull_pool", "hull_top128")):
        va = G(a) if a in ("best1_top128", "bestm128", "hull_top128") else (
            Mc[:, :M].min(1) if a == "bestm_top75" else (Mc.min(1) if a == "bestm_pool" else C(a)))
        vb = G(b) if b in ("best1_top128", "bestm128", "hull_top128") else (
            Mc[:, :M].min(1) if b == "bestm_top75" else (Mc.min(1) if b == "bestm_pool" else C(b)))
        out.setdefault("prefix_contrasts", {})[a + " vs " + b] = ST.compare(
            va, vb, folds, names=pdbs, label="ORACLE %s vs ORACLE %s (POINT CLOUD)" % (a, b))
    ST.save_atomic(P128_JSON, out, module_file=__file__)
    print("== rung 9 (the top-128 prefix: the deployed quantum stage's ENTIRE field of view). "
          "The prefix is native-free; EVERY RMSD below is ORACLE.")
    print("   production %.4f; the top-75 is %.0f%% inside the top-128 prefix by construction; "
          "ORACLE best prefix-m inside 128: median m %.0f (10/25/50/75/90 %s), at m=128 on %.0f%% of targets; "
          "hull support median %.0f" % (out["prod_mean"], 100 * out["overlap_mean"], out["m128_median"],
                                        ["%.0f" % q for q in out["m128_quartiles"]],
                                        100 * out["frac_m128_at_128"], out["hull_support_median"]))
    print("   | class | top-75 | TOP-128 | K=500 |")
    for cls, keys in (("best single member", ("best1_top75", "best1_top128", "best1_pool")),
                      ("best prefix-m average", ("bestm_top75", "bestm128", "bestm_pool")),
                      ("convex hull", ("hull_top75", "hull_top128", "hull_pool"))):
        print("   | %s | %.4f | %.4f | %.4f |" % (cls, trip[keys[0]]["mean"], trip[keys[1]]["mean"],
                                                  trip[keys[2]]["mean"]))
    for k, v in out["prefix_contrasts"].items():
        print(ST.fmt(v))
    print("p128:", P128_JSON)
    return out


# ================================================================== per-target point cloud
def cloud_row(pdb, store=True, verbose=True):
    t0 = time.time()
    cand, dis, top, order, dg = load_pool(pdb)
    W = cand.W; n = cand.n; nat = cand.nat_ca; rr = cand.oracle_rr           # ORACLE labels
    u = I.load_univ(pdb)
    fail = pdb in I.FAIL18
    P = I.pairwise_rmsd(W)
    structs: Dict[str, np.ndarray] = {}
    row: Dict = dict(pdb=pdb, n=n, fold=int(cand.fold), fail18=bool(fail), items={})

    # production (the anchor and rung 6's C)
    C, b = avg_of(W, top, P)
    C_ref, _ = I.coordinate_average(W[top])
    assert np.abs(C - C_ref).max() < 1e-9, "avg_of != coordinate_average"
    f28 = os.path.join(ROOT, "s27", "results", "s28_A_structs", f"{pdb}.npz")
    if os.path.exists(f28):
        with np.load(f28) as z:
            row["prod_dev_vs_s28"] = float(np.abs(np.asarray(z["prod"], float) - C).max())
    row["items"]["prod"] = dict(rmsd_cloud=float(I.ca_rmsd(C, nat)), **struct_diag(C))     # ORACLE
    structs["prod"] = C

    # rung 1 -- best single member                                                     # ORACLE
    for name, idx in (("best1_top75", np.sort(top)), ("best1_pool", np.arange(cand.k))):
        r = rr[idx]; m = r.min(); tie = np.isclose(r, m, atol=1e-9)
        j = int(idx[np.flatnonzero(tie)[0]])
        X = W[j]
        row["items"][name] = dict(rmsd_cloud=float(I.ca_rmsd(X, nat)), rr=float(m), n_tied=int(tie.sum()),
                                  member=j, **struct_diag(X))
        structs[name] = X

    # rung 2 -- the m-ladder of the DIS order                                          # ORACLE
    curve = oracle_m_curve(W, order, P, nat)
    assert abs(curve[M - 1] - row["items"]["prod"]["rmsd_cloud"]) < 1e-9, "m=75 != production"
    m_best = int(np.argmin(curve)) + 1
    Xm, _ = avg_of(W, order[:m_best], P)
    row["m_curve"] = curve.tolist()
    row["items"]["bestm"] = dict(rmsd_cloud=float(curve[m_best - 1]), m=m_best, **struct_diag(Xm))
    structs["bestm"] = Xm

    # rung 3 -- best basin average                                                     # ORACLE
    sets = {"top75": np.sort(top), "pool": np.arange(cand.k)}
    for S, idx in sets.items():
        for k in BASIN_K:
            r, X, size, n_found, n_elig = oracle_best_basin(W, idx, P, nat, k)
            name = "basin_%s_k%d" % (S, k)
            row["items"][name] = dict(rmsd_cloud=r, size=size, n_basins=n_found, n_eligible=n_elig, **struct_diag(X))
            structs[name] = X

    # rungs 4 and 5 -- sparse convex and the hull, in the S28 frame                    # ORACLE
    from s27 import s28_A_amp as A
    frame = A.Frame(W, np.sort(top))
    natp0 = superpose_one(nat, frame.ref)
    for S, idx in sets.items():
        Wf = frame.Wf[idx]
        r, X, w = oracle_hull(Wf, nat, n, natp0=natp0)
        row["items"]["hull_%s" % S] = dict(rmsd_cloud=r, n_support=int((w > 1e-9).sum()), **struct_diag(X))
        structs["hull_%s" % S] = X
        sp = oracle_sparse_greedy(Wf, nat, n, natp0=natp0)
        for s, (r, X, sup, w) in sp.items():
            name = "sparse_%s_s%d" % (S, s)
            row["items"][name] = dict(rmsd_cloud=r, support=[int(idx[q]) for q in sup], w=[float(x) for x in w],
                                      **struct_diag(X))
            structs[name] = X

    # rung 6 -- the typicality axis, both blind definitions, gated against S24's artefacts
    row["axis"] = {}
    for which in BLIND_DEFS:
        B, mem, meta = blind_cloud(pdb, u, dg, which)
        got, want, CA24 = s24_gate_values(pdb, u, dg, B, which)
        dev = {k: abs(got[k] - want[k]) for k in want}
        if max(dev.values()) > 1e-6:
            raise RuntimeError("S24 gate failed for %s on %s: got %s want %s" % (which, pdb, got, want))
        assert np.abs(CA24 - C).max() < 1e-6, "S24's A average != production cloud"
        rec, uvec = axis_probe(C, B, nat, pdb, which)
        rec.update(gate=dict(got=got, want=want, max_dev=max(dev.values())), members=[int(x) for x in mem], **meta)
        row["axis"][which] = rec
        structs["blind_%s" % which] = B
        structs["u_%s" % which] = uvec.reshape(n, 3)

    row["secs_cloud"] = time.time() - t0
    if store:
        save_structs(pdb, structs, cand.seq, cand.fold, cand.n)
    if verbose:
        it = row["items"]
        print("  %s n=%d prod %.3f | best1 %.3f/%.3f | bestm %.3f (m=%d) | basin75 %.3f basin500 %.3f | "
              "sparse75 s2 %.3f s20 %.3f | sparse500 s2 %.3f s20 %.3f | hull %.3f/%.3f | cos LIB75 %+.3f BPRIME %+.3f | %.1fs"
              % (pdb, n, it["prod"]["rmsd_cloud"], it["best1_top75"]["rmsd_cloud"], it["best1_pool"]["rmsd_cloud"],
                 it["bestm"]["rmsd_cloud"], it["bestm"]["m"],
                 min(it["basin_top75_k%d" % k]["rmsd_cloud"] for k in BASIN_K),
                 min(it["basin_pool_k%d" % k]["rmsd_cloud"] for k in BASIN_K),
                 it["sparse_top75_s2"]["rmsd_cloud"], it["sparse_top75_s20"]["rmsd_cloud"],
                 it["sparse_pool_s2"]["rmsd_cloud"], it["sparse_pool_s20"]["rmsd_cloud"],
                 it["hull_top75"]["rmsd_cloud"], it["hull_pool"]["rmsd_cloud"],
                 row["axis"]["LIB75"]["cos_uv"], row["axis"]["BPRIME"]["cos_uv"], row["secs_cloud"]), flush=True)
    return row, structs


# ============================================================================ storage
def replace_retry(tmp, f, tries=8, wait=0.5):
    for k in range(int(tries)):
        try:
            os.replace(tmp, f); return
        except PermissionError:
            if k == tries - 1:
                raise
            time.sleep(wait)


def save_structs(pdb, structs, seq, fold, n, merge=True):
    os.makedirs(STRUCTS, exist_ok=True)
    f = os.path.join(STRUCTS, f"{pdb}.npz")
    old = {}
    if merge and os.path.exists(f):
        with np.load(f) as z:
            old = {k: np.array(z[k]) for k in z.files}
    old.update({k: np.asarray(v, float) for k, v in structs.items()})
    old["seq"] = np.array(str(seq)); old["fold"] = np.array(int(fold)); old["n"] = np.array(int(n))
    tmp = f + f".tmp{os.getpid()}.npz"
    np.savez_compressed(tmp, **old)
    replace_retry(tmp, f)


def load_structs(pdb):
    f = os.path.join(STRUCTS, f"{pdb}.npz")
    with np.load(f) as z:
        return {k: np.array(z[k]) for k in z.files}


def jsonl_rows(path, key=("pdb",)):
    done = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except Exception:
                    continue
                done[tuple(r[k] for k in key)] = r
    return done


def append_row(path, row):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)) + "\n")


def all_pdbs():
    return [t["pdb"] for t in I.targets()]


# ============================================================================= phases
def phase_cloud(pdbs):
    done = jsonl_rows(CLOUD_ROWS)
    t0 = time.time()
    for k, pdb in enumerate(pdbs):
        if (pdb,) in done:
            continue
        row, _ = cloud_row(pdb)
        append_row(CLOUD_ROWS, row)
        print("  [cloud %d/%d] elapsed %.1f min" % (k + 1, len(pdbs), (time.time() - t0) / 60), flush=True)
    print("cloud rows:", CLOUD_ROWS)


def lfo_choices(curves, folds, t_grid=T_GRID, only_folds=None):
    """Leave-fold-out step: for each fold f, t_f = argmin over the grid of the mean curve over the
    targets NOT in f.  `curves` (n_targets, n_t) are ORACLE per-target curves; the choice for a
    target never reads its own fold, which is what `only_folds` lets the NaN-poison test show.
    A non-finite TRAINING slice is refused.  Returns (t per target, t per fold, ties per fold)."""
    curves = np.asarray(curves, float); folds = np.asarray(folds)
    want = sorted(set(folds.tolist())) if only_folds is None else [int(f) for f in only_folds]
    t_of = np.full(len(folds), np.nan); t_fold = {}; ties = {}
    for f in want:
        sl = curves[folds != f]
        if not np.isfinite(sl).all():
            raise ValueError("non-finite training curves for fold %s" % f)
        train = sl.mean(0)
        m = train.min(); tie = np.isclose(train, m, atol=1e-12)
        j = int(np.flatnonzero(tie)[0])
        t_fold[int(f)] = float(t_grid[j]); ties[int(f)] = int(tie.sum())
        t_of[folds == f] = t_grid[j]
    return t_of, t_fold, ties


def rand18_null(d, fail_mask, n_draw=20000, seed_parts=("s29O", "rand18")):
    """The FAIL18 mean of a per-target difference against 20,000 random 18-subsets of the 126
    (without replacement); one-sided p in the observed direction (S28-L40's construction)."""
    d = np.asarray(d, float); fm = np.asarray(fail_mask, bool)
    k = int(fm.sum()); n = len(d)
    obs = float(d[fm].mean()); other = float(d[~fm].mean())
    rng = SD.stable_rng(*seed_parts, salt=SALT)
    draws = np.array([d[rng.choice(n, k, replace=False)].mean() for _ in range(int(n_draw))])
    direction = 1.0 if obs >= draws.mean() else -1.0
    p = float(((direction * draws) >= (direction * obs)).mean())
    return dict(fail18_mean=obs, other_mean=other, k=k, n_draw=int(n_draw), p_one_sided=p,
                null_p2_5=float(np.percentile(draws, 2.5)), null_p97_5=float(np.percentile(draws, 97.5)),
                null_mean=float(draws.mean()))


def strata(a, b, pdbs, folds, label):
    """ST.compare on the 126, on FAIL18 and on the 108, plus the random-18 null on the difference."""
    a = np.asarray(a, float); b = np.asarray(b, float); folds = np.asarray(folds)
    fm = np.array([p in I.FAIL18 for p in pdbs])
    out = dict(all=ST.compare(a, b, folds, names=pdbs, label=label))
    out["fail18"] = ST.compare(a[fm], b[fm], folds[fm], names=[p for p, f in zip(pdbs, fm) if f], label=label + " [FAIL18]")
    out["other108"] = ST.compare(a[~fm], b[~fm], folds[~fm], names=[p for p, f in zip(pdbs, fm) if not f], label=label + " [108]")
    out["rand18"] = rand18_null(a - b, fm)
    return out


def phase_lfo():
    rows = jsonl_rows(CLOUD_ROWS)
    pdbs = all_pdbs()
    missing = [p for p in pdbs if (p,) not in rows]
    if missing:
        raise RuntimeError("lfo needs all 126 cloud rows; missing %d: %s" % (len(missing), missing[:5]))
    folds = ST.pinned_folds(pdbs)
    prod = np.array([rows[(p,)]["items"]["prod"]["rmsd_cloud"] for p in pdbs])
    out = dict(t_grid=[float(t) for t in T_GRID], defs={}, provenance=ST.provenance(__file__))
    for which in BLIND_DEFS:
        R = np.array([rows[(p,)]["axis"][which]["curve"] for p in pdbs])         # ORACLE curves
        t_of, t_fold, ties = lfo_choices(R, folds)
        mean_curve = R.mean(0)
        j_glob = int(np.argmin(mean_curve))
        arm = np.array([R[i, int(np.argmin(np.abs(T_GRID - t_of[i])))] for i in range(len(pdbs))])
        assert np.allclose(R[:, int(np.argmin(np.abs(T_GRID)))], prod), "t=0 column != production"
        cos = np.array([rows[(p,)]["axis"][which]["cos_uv"] for p in pdbs])
        cr = np.array([rows[(p,)]["axis"][which]["cos_rand_mean"] for p in pdbs])
        cra = np.array([rows[(p,)]["axis"][which]["cos_rand_absmean"] for p in pdbs])
        d = dict(
            t_fold=t_fold, ties=ties, t_of_target={p: float(t) for p, t in zip(pdbs, t_of)},
            mean_curve=mean_curve.tolist(),
            oracle_global_t=float(T_GRID[j_glob]), oracle_global_mean=float(mean_curve[j_glob]),
            oracle_per_target_mean=float(R.min(1).mean()),
            oracle_per_target_bok=ST.best_of_k_within(R, seed_parts=("s29O", which, "tgrid")),
            cos=dict(mean=float(cos.mean()), median=float(np.median(cos)), se=float(cos.std(ddof=1) / math.sqrt(len(cos))),
                     vs_rand_signed=ST.compare(cos, cr, folds, names=pdbs, label="cos(u,v) vs random-field signed cos [%s]" % which),
                     vs_rand_abs=ST.compare(cos, cra, folds, names=pdbs, label="cos(u,v) vs random-field |cos| [%s]" % which),
                     fail18_mean=float(cos[[p in I.FAIL18 for p in pdbs]].mean()),
                     other108_mean=float(cos[[p not in I.FAIL18 for p in pdbs]].mean()),
                     rand18=rand18_null(cos, [p in I.FAIL18 for p in pdbs], seed_parts=("s29O", "rand18cos", which)),
                     n_positive=int((cos > 0).sum()), n_above_rand_abs=int((cos > cra).sum())),
            lfo_cloud=strata(arm, prod, pdbs, folds, "lfo_%s vs prod (POINT CLOUD, DEPLOYABLE step)" % which),
            blind_mean=float(np.mean([rows[(p,)]["axis"][which]["rmsd_blind"] for p in pdbs])),
            arm_mean=float(arm.mean()), prod_mean=float(prod.mean()),
            arm_cloud={p: float(v) for p, v in zip(pdbs, arm)},
            cos_rand_analytic_absmean=float(np.mean([rows[(p,)]["axis"][which]["cos_rand_analytic_absmean"] for p in pdbs])),
        )
        out["defs"][which] = d
        # store the deployable structures
        for i, p in enumerate(pdbs):
            z = load_structs(p)
            X = lfo_structure(z["prod"], z["u_%s" % which], t_of[i])
            save_structs(p, {"lfo_%s" % which: X}, str(z["seq"]), int(z["fold"]), int(z["n"]))
        print("== %s: t per fold %s (ties %s); ORACLE global t %.1f -> %.4f; ORACLE per-target %.4f; LFO arm %.4f vs prod %.4f"
              % (which, t_fold, ties, d["oracle_global_t"], d["oracle_global_mean"], d["oracle_per_target_mean"],
                 d["arm_mean"], d["prod_mean"]), flush=True)
        print("   cos(u,v) mean %+.4f median %+.4f ; random signed %+.4f ; random |cos| %.4f ; FAIL18 %+.4f / 108 %+.4f"
              % (d["cos"]["mean"], d["cos"]["median"], cr.mean(), cra.mean(), d["cos"]["fail18_mean"], d["cos"]["other108_mean"]))
        print(ST.fmt(d["lfo_cloud"]["all"]))
    ST.save_atomic(LFO_JSON, out, module_file=__file__)
    print("lfo:", LFO_JSON)
    return out


def chain_item(pdb, item):
    z = load_structs(pdb)
    if item not in z:
        return dict(pdb=pdb, item=item, rmsd_chain=float("nan"), missing=True)
    C = np.asarray(z[item], float)
    u = I.load_univ(pdb); nat = u["nat_ca"]                                    # ORACLE scoring
    if not np.isfinite(C).all():
        return dict(pdb=pdb, item=item, rmsd_chain=float("nan"), undefined=True)
    t1 = time.time()
    pr = I.project(C, str(z["seq"]), int(z["fold"]))
    ca = np.asarray(pr["ca"], float)
    return dict(pdb=pdb, item=item, rmsd_chain=float(I.ca_rmsd(ca, nat)), rmsd_cloud=float(I.ca_rmsd(C, nat)),
                secs=time.time() - t1, **struct_diag(ca))


def chain_rows_path(shard=None):
    return CHAIN_ROWS if shard is None else CHAIN_ROWS.replace(".jsonl", "_shard%d.jsonl" % int(shard))


def all_chain_rows():
    """Every chain row from the unsharded file and from every shard file."""
    import glob as _glob
    done = {}
    for f in sorted(_glob.glob(CHAIN_ROWS.replace(".jsonl", "*.jsonl"))):
        done.update(jsonl_rows(f, key=("pdb", "item")))
    return done


def phase_chain(pdbs, groups="ABCD", shard=None, n_shards=1):
    #: every shard READS every rows file (so nothing is computed twice) and APPENDS only to its own
    #: (so no two processes ever write the same file). Work splits by target index: a kill costs one
    #: item, a resume skips whatever any shard already banked.
    done = all_chain_rows()
    path = chain_rows_path(shard)
    if shard is not None:
        pdbs = [p for k, p in enumerate(pdbs) if k % int(n_shards) == int(shard)]
    work = []
    for g in groups:
        if g == "A":
            work += [(p, it) for it in CHAIN_GROUPS["A"] for p in pdbs]           # rung 6 first, all targets
        else:
            work += [(p, it) for p in pdbs for it in CHAIN_GROUPS[g]]
    todo = [w for w in work if w not in done]
    print("chain shard %s/%s: %d targets, %d items, %d done, %d to do -> %s"
          % (shard, n_shards, len(pdbs), len(work), len(work) - len(todo), len(todo),
             os.path.basename(path)), flush=True)
    t0 = time.time()
    for k, (pdb, item) in enumerate(todo):
        r = chain_item(pdb, item)
        append_row(path, r)
        if (k + 1) % 20 == 0 or k == len(todo) - 1:
            print("  [chain %d/%d] %s %s chain %.3f cloud %.3f (%.1fs)  elapsed %.1f min"
                  % (k + 1, len(todo), pdb, item, r.get("rmsd_chain", float("nan")), r.get("rmsd_cloud", float("nan")),
                     r.get("secs", 0.0), (time.time() - t0) / 60), flush=True)
    print("chain rows:", path)


# ============================================================================ analysis
def _fold_ci_mean(x, folds, n_boot=4000, seed_parts=("s29O", "foldci")):
    x = np.asarray(x, float); folds = np.asarray(folds)
    F = np.array(sorted(set(folds.tolist())))
    rng = SD.stable_rng(*seed_parts, salt=SALT)
    b = np.array([np.concatenate([x[folds == q] for q in rng.choice(F, len(F), replace=True)]).mean() for _ in range(n_boot)])
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def analyse(write=True):
    rows = jsonl_rows(CLOUD_ROWS)
    chain = all_chain_rows()
    pdbs = [p for p in all_pdbs() if (p,) in rows]
    folds = ST.pinned_folds(pdbs)
    fm = np.array([p in I.FAIL18 for p in pdbs])
    out = dict(n_cloud=len(pdbs), n_chain_items=len(chain), items={}, provenance=ST.provenance(__file__))
    prod_chain = np.array([chain.get((p, "prod"), {}).get("rmsd_chain", np.nan) for p in pdbs])
    lfo = json.load(open(LFO_JSON)) if os.path.exists(LFO_JSON) else None
    for it in ALL_ITEMS:
        rc = np.full(len(pdbs), np.nan)
        if not it.startswith("lfo_"):
            rc = np.array([rows[(p,)]["items"].get(it, {}).get("rmsd_cloud", np.nan) for p in pdbs])
            if it in ("best1_top128", "bestm128", "hull_top128") and os.path.exists(P128_ROWS):
                p128 = jsonl_rows(P128_ROWS)
                rc = np.array([p128.get((p,), {}).get(it, np.nan) for p in pdbs])
        rh = np.array([chain.get((p, it), {}).get("rmsd_chain", np.nan) for p in pdbs])
        rhc = np.array([chain.get((p, it), {}).get("rmsd_cloud", np.nan) for p in pdbs])
        if it.startswith("lfo_") and lfo is not None:
            rc = np.array([lfo["defs"][it[4:]]["arm_cloud"].get(p, np.nan) for p in pdbs])
        ok_c = np.isfinite(rc); ok_h = np.isfinite(rh)
        d = dict(cloud=dict(n=int(ok_c.sum()), mean=float(np.nanmean(rc)) if ok_c.any() else None,
                            median=float(np.nanmedian(rc)) if ok_c.any() else None,
                            fail18=float(np.nanmean(rc[fm])) if ok_c[fm].any() else None,
                            other108=float(np.nanmean(rc[~fm])) if ok_c[~fm].any() else None,
                            frac_under2=float(np.nanmean(rc < 2.0)) if ok_c.any() else None),
                 chain=dict(n=int(ok_h.sum()), mean=float(np.nanmean(rh)) if ok_h.any() else None,
                            median=float(np.nanmedian(rh)) if ok_h.any() else None,
                            fail18=float(np.nanmean(rh[fm])) if ok_h[fm].any() else None,
                            other108=float(np.nanmean(rh[~fm])) if ok_h[~fm].any() else None,
                            frac_under2=float(np.nanmean(rh[ok_h] < 2.0)) if ok_h.any() else None))
        if ok_h.sum() >= 5:
            pr = (rh - rhc)[ok_h]
            d["price"] = dict(mean=float(pr.mean()), fold_ci=_fold_ci_mean(pr, folds[ok_h], seed_parts=("s29O", "price", it)),
                              fail18=float(pr[fm[ok_h]].mean()) if fm[ok_h].any() else None,
                              other108=float(pr[~fm[ok_h]].mean()) if (~fm[ok_h]).any() else None)
            if it != "prod" and np.isfinite(prod_chain[ok_h]).all():
                d["vs_prod_chain"] = ST.compare(rh[ok_h], prod_chain[ok_h], folds[ok_h], names=[p for p, o in zip(pdbs, ok_h) if o],
                                                label="%s vs prod (BUILT CHAIN)%s" % (it, "" if it.startswith("lfo_") else " [ORACLE]"))
        out["items"][it] = d
    # rung 2: the order statistic over m
    Mc = np.array([rows[(p,)]["m_curve"] for p in pdbs])
    out["rung2"] = dict(mean_curve_at=[1, 5, 10, 25, 50, 75, 100, 150, 200, 300, 500],
                        mean_curve=[float(Mc[:, m - 1].mean()) for m in [1, 5, 10, 25, 50, 75, 100, 150, 200, 300, 500]],
                        global_best_m=int(np.argmin(Mc.mean(0))) + 1, global_best_mean=float(Mc.mean(0).min()),
                        per_target_best_mean=float(Mc.min(1).mean()),
                        m_best_median=float(np.median([rows[(p,)]["items"]["bestm"]["m"] for p in pdbs])),
                        m_best_quartiles=[float(q) for q in np.percentile([rows[(p,)]["items"]["bestm"]["m"] for p in pdbs], [10, 25, 50, 75, 90])],
                        bok=ST.best_of_k_within(Mc, seed_parts=("s29O", "mcurve")))
    # rung 3: which cell wins, basin sizes
    cells = ["basin_%s_k%d" % (S, k) for S in SETS for k in BASIN_K]
    Bc = np.array([[rows[(p,)]["items"][c]["rmsd_cloud"] for c in cells] for p in pdbs])
    win = np.nanargmin(Bc, 1)
    out["rung3"] = dict(cells=cells, cell_means=[float(np.nanmean(Bc[:, j])) for j in range(len(cells))],
                        best_over_cells_mean=float(np.nanmin(Bc, 1).mean()),
                        win_counts={cells[j]: int((win == j).sum()) for j in range(len(cells))},
                        win_size_median=float(np.median([rows[(p,)]["items"][cells[w]]["size"] for p, w in zip(pdbs, win)])),
                        size_by_cell={c: float(np.median([rows[(p,)]["items"][c]["size"] for p in pdbs])) for c in cells},
                        bok=ST.best_of_k_within(np.nan_to_num(Bc, nan=99.0), seed_parts=("s29O", "basins")))
    # rung 6
    if os.path.exists(LFO_JSON):
        out["lfo"] = json.load(open(LFO_JSON))
        for which in BLIND_DEFS:
            it = "lfo_%s" % which
            rh = np.array([chain.get((p, it), {}).get("rmsd_chain", np.nan) for p in pdbs])
            ok = np.isfinite(rh) & np.isfinite(prod_chain)
            if ok.sum() >= 5:
                out["lfo"]["defs"][which]["lfo_chain"] = strata(rh[ok], prod_chain[ok], [p for p, o in zip(pdbs, ok) if o], folds[ok],
                                                                "lfo_%s vs prod (BUILT CHAIN, DEPLOYABLE step)" % which)
                out["lfo"]["defs"][which]["lfo_chain_n"] = int(ok.sum())
    if os.path.exists(PC1_JSON):
        out["pc1"] = json.load(open(PC1_JSON))
    if os.path.exists(P128_JSON):
        out["p128"] = json.load(open(P128_JSON))
    out["text"] = render_table(out)
    if write:
        ST.save_atomic(SUMMARY_JSON, out, module_file=__file__)
    print(out["text"])
    return out


def _f(x, w=6):
    return ("%%%d.3f" % w) % x if x is not None and np.isfinite(x) else " " * (w - 1) + "-"


def render_table(out):
    L = []
    L.append("| item (all ORACLE except lfo_*) | cloud mean | cloud FAIL18 | cloud 108 | cloud <2A | chain n | chain mean | chain FAIL18 | chain 108 | chain <2A | price (chain-cloud) [fold CI] |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for it in ALL_ITEMS:
        d = out["items"][it]; c = d["cloud"]; h = d["chain"]; pr = d.get("price")
        L.append("| %s | %s | %s | %s | %s | %d | %s | %s | %s | %s | %s |" % (
            it, _f(c["mean"]), _f(c["fail18"]), _f(c["other108"]), _f(c["frac_under2"]), h["n"], _f(h["mean"]), _f(h["fail18"]),
            _f(h["other108"]), _f(h["frac_under2"]),
            ("%+.3f [%+.3f, %+.3f]" % (pr["mean"], pr["fold_ci"][0], pr["fold_ci"][1])) if pr else "-"))
    r2 = out["rung2"]
    L.append("")
    L.append("rung 2 (ORACLE m-ladder of the DIS order): mean curve at m=%s -> %s; ORACLE global best m=%d (%.4f); ORACLE per-target best m: mean %.4f, median m %.0f (10/25/50/75/90 pct %s); best-of-500 pricing: observed gain %+.4f, across-target null %+.4f (%.0f%% accounted), split-half transfer %+.4f (%.0f%% of the oracle), k_eff %.1f"
             % (r2["mean_curve_at"], ["%.3f" % v for v in r2["mean_curve"]], r2["global_best_m"], r2["global_best_mean"],
                r2["per_target_best_mean"], r2["m_best_median"], ["%.0f" % q for q in r2["m_best_quartiles"]],
                r2["bok"]["observed_gain"], r2["bok"]["null_across_targets"], 100 * r2["bok"]["share_accounted"],
                r2["bok"]["split_half"], 100 * r2["bok"]["split_half_frac"], r2["bok"]["k_eff"]))
    r3 = out["rung3"]
    L.append("rung 3 (ORACLE best basin): cell means %s; best over the 10 cells %.4f (an order statistic: split-half transfer %+.4f, %.0f%% of the oracle); win counts %s; winning basin size median %.0f; size medians by cell %s"
             % ({c: "%.3f" % v for c, v in zip(r3["cells"], r3["cell_means"])}, r3["best_over_cells_mean"], r3["bok"]["split_half"],
                100 * r3["bok"]["split_half_frac"], r3["win_counts"], r3["win_size_median"], {c: "%.0f" % v for c, v in r3["size_by_cell"].items()}))
    if "lfo" in out:
        for which in BLIND_DEFS:
            d = out["lfo"]["defs"][which]
            L.append("")
            L.append("rung 6 [%s]: cos(u,v) ORACLE mean %+.4f (SE %.4f, median %+.4f), random-field signed %+.4f, random |cos| %.4f (analytic %.4f); FAIL18 %+.4f / 108 %+.4f (random-18 p %.3f); positive on %d/126, above random |cos| on %d/126"
                     % (which, d["cos"]["mean"], d["cos"]["se"], d["cos"]["median"], d["cos"]["vs_rand_signed"]["mean_b"],
                        d["cos"]["vs_rand_abs"]["mean_b"], d["cos_rand_analytic_absmean"],
                        d["cos"]["fail18_mean"], d["cos"]["other108_mean"], d["cos"]["rand18"]["p_one_sided"], d["cos"]["n_positive"], d["cos"]["n_above_rand_abs"]))
            L.append("   blind average ORACLE mean %.4f; ORACLE global step t*=%.1f -> %.4f; ORACLE per-target step mean %.4f (best-of-31: split-half %+.4f, %.0f%%); LEAVE-FOLD-OUT t per fold %s -> arm %.4f vs prod %.4f (point cloud)"
                     % (d["blind_mean"], d["oracle_global_t"], d["oracle_global_mean"], d["oracle_per_target_mean"],
                        d["oracle_per_target_bok"]["split_half"], 100 * d["oracle_per_target_bok"]["split_half_frac"], d["t_fold"], d["arm_mean"], d["prod_mean"]))
            L.append(ST.fmt(d["cos"]["vs_rand_signed"]))
            L.append(ST.fmt(d["lfo_cloud"]["all"]))
            if "lfo_chain" in d:
                L.append(ST.fmt(d["lfo_chain"]["all"]))
                L.append("   chain strata: FAIL18 %+.4f (%.2fx MDE, %d/%d folds) / 108 %+.4f (%.2fx MDE); random-18 null p %.3f [null 2.5/97.5 %+.4f, %+.4f]"
                         % (d["lfo_chain"]["fail18"]["effect"], d["lfo_chain"]["fail18"]["effect_over_mde"], d["lfo_chain"]["fail18"]["folds_same_sign"],
                            d["lfo_chain"]["fail18"]["n_folds"], d["lfo_chain"]["other108"]["effect"], d["lfo_chain"]["other108"]["effect_over_mde"],
                            d["lfo_chain"]["rand18"]["p_one_sided"], d["lfo_chain"]["rand18"]["null_p2_5"], d["lfo_chain"]["rand18"]["null_p97_5"]))
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["cloud", "lfo", "chain", "pc1", "p128", "analyse"])
    ap.add_argument("--pdbs", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--groups", default="ABCD")
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--n-shards", type=int, default=1)
    a = ap.parse_args()
    pdbs = all_pdbs()
    if a.pdbs:
        pdbs = [p for p in a.pdbs.split(",") if p]
    if a.limit:
        pdbs = pdbs[:a.limit]
    os.makedirs(RESULTS, exist_ok=True)
    if a.mode == "cloud":
        phase_cloud(pdbs)
    elif a.mode == "lfo":
        phase_lfo()
    elif a.mode == "chain":
        phase_chain(pdbs, groups=a.groups, shard=a.shard, n_shards=a.n_shards)
    elif a.mode == "pc1":
        phase_pc1(pdbs)
    elif a.mode == "p128":
        phase_p128(pdbs)
    else:
        analyse()


if __name__ == "__main__":
    main()
