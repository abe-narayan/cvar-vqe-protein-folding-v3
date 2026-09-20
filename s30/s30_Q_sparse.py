#!/usr/bin/env python
"""s30/s30_Q_sparse.py -- S30 lane Q, measurement 1 of `s30/PREREG_S30_Q_sparse.md`.

The sparse weighted readout's ceiling as a function of (a) support size s, (b) HOW the support is
chosen, and (c) the BIT PRICE of support and weights counted in the same currency.

EVERY ARM HERE IS ORACLE except the SCORE-PREFIX-UNIFORM and DIVERSITY-UNIFORM arms, which read no
native. ORACLE reads happen only inside functions whose names begin `oracle_` or under `# ORACLE`
markers. No deployable parameter is chosen anywhere in this file.

Basis: POINT CLOUD. The built chain is the endpoint (S30 contract rule 1) and is run separately on
the decisive arms; the ORACLE/ORACLE corner of this family is ALREADY chained by S29 lane O.

    python s30/s30_Q_sparse.py cloud [--limit N] [--pdbs A,B]
    python s30/s30_Q_sparse.py analyse
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time
from typing import Dict, List, Sequence, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                    # noqa: E402
from s24 import stats_lib as ST                    # noqa: E402
from s29 import s29_O_ladder as O                  # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s30_Q_sparse_rows.jsonl")
SUMMARY = os.path.join(RESULTS, "s30_Q_sparse.json")
O_CLOUD = os.path.join(ROOT, "s29", "results", "s29_O_cloud_rows.jsonl")

SALT = "s30Q"
K = 500
T_LADDER = (4, 8, 16, 32, 64, 128, 256, 500)       # score-order prefixes the support is drawn from
S_LADDER = (2, 3, 5)                               # support sizes on the T ladder
S_FULL = (2, 3, 5, 10, 20)                         # support sizes on the factorial (sets top75/pool)
L_GRID = (1, 2, 3, 4, 8, 16, 32)                   # weight-grid resolutions (multiples of 1/L)
N_RAND = 8                                         # random-support draws (MEAN is the arm)
B_BITS = tuple(range(0, 10))                       # reference argmin over top-2**B
EXHAUST_T = 64                                     # exhaustive s=2 check bounding the greedy gap
ENUM_MAX = 4000                                    # enumerate the weight grid below this many points


# ------------------------------------------------------------------ bit ledger
def log2C(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("nan")
    return float(math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2.0)


def support_bits(pool_size: int, s: int) -> float:
    """Bits to name WHICH s of `pool_size`."""
    return log2C(int(pool_size), int(s))


def weight_bits(s: int, L) -> float:
    """Bits to name a weight vector on the simplex grid with denominator L (lattice-point count).
    L None/inf -> continuous, unbounded; L == 0 -> uniform, free."""
    if L is None:
        return float("inf")
    if L == 0:
        return 0.0
    return log2C(int(L) + int(s) - 1, int(s) - 1)


# ------------------------------------------------------------------ geometry helpers
def rmsd_batch(X, nat):
    """CA-RMSD of each X[k] (n,3) to nat, Kabsch-aligned."""
    return np.asarray(I.kabsch_rmsd_batch(np.asarray(X, float), np.asarray(nat, float)), float)


def combos_on_simplex(s: int, L: int) -> np.ndarray:
    """All weight vectors that are multiples of 1/L summing to 1 (lattice points), shape (G, s)."""
    out = []
    for cut in itertools.combinations(range(L + s - 1), s - 1):
        prev, row = -1, []
        for c in cut:
            row.append(c - prev - 1)
            prev = c
        row.append(L + s - 2 - prev)
        out.append(row)
    return np.asarray(out, float) / float(L)


def oracle_grid_weights(Wf_sup, nat, n, L, w0=None):
    """ORACLE. Best weight vector on the 1/L simplex lattice over a GIVEN support.
    Exhaustive when the lattice is small; otherwise largest-remainder rounding of `w0` followed by
    a pairwise 1/L local search. Returns (rmsd, w)."""
    s = len(Wf_sup)
    G = log2C(L + s - 1, s - 1)
    if G <= math.log2(ENUM_MAX):
        Wg = combos_on_simplex(s, L)                                    # (G, s)
        X = (Wg @ Wf_sup).reshape(len(Wg), n, 3)
        r = rmsd_batch(X, nat)                                          # ORACLE
        j = int(np.argmin(r))
        return float(r[j]), Wg[j].copy()
    # largest-remainder rounding of the continuous optimum, then local search
    w0 = np.full(s, 1.0 / s) if w0 is None else np.asarray(w0, float)
    q = np.floor(w0 * L).astype(int)
    rem = L - q.sum()
    if rem > 0:
        q[np.argsort(-(w0 * L - q))[:rem]] += 1
    w = q.astype(float) / L
    best = float(rmsd_batch((w @ Wf_sup).reshape(1, n, 3), nat)[0])     # ORACLE
    for _ in range(60):
        cand_w, moves = [], []
        for a in range(s):
            if q[a] <= 0:
                continue
            for b in range(s):
                if a == b:
                    continue
                q2 = q.copy(); q2[a] -= 1; q2[b] += 1
                cand_w.append(q2 / float(L)); moves.append(q2)
        if not cand_w:
            break
        Wc = np.asarray(cand_w, float)
        r = rmsd_batch((Wc @ Wf_sup).reshape(len(Wc), n, 3), nat)       # ORACLE
        j = int(np.argmin(r))
        if r[j] < best - 1e-12:
            best = float(r[j]); q = moves[j]; w = Wc[j]
        else:
            break
    return best, np.asarray(w, float)


def oracle_greedy_support(Wf, idx, nat, n, s_max, natp0, rounds=O.SPARSE_ROUNDS):
    """ORACLE. Lane O's greedy forward support selection restricted to candidate rows `idx`.
    Returns {s: (rmsd, support_global, w)} for every s = 1..s_max (nested path)."""
    idx = np.asarray(idx, int)
    A = Wf[idx]
    k = len(A)
    natp = np.asarray(natp0, float).ravel()
    support: List[int] = []
    out = {}
    w = None
    for step in range(1, min(int(s_max), k) + 1):
        best = None
        for j in range(k):
            if j in support:
                continue
            cols = support + [j]
            wj = O.convex_nnls(A[cols].T, natp)                          # ORACLE (natp is the native)
            res = float(np.linalg.norm(wj @ A[cols] - natp))
            if best is None or res < best[0] - 1e-12:
                best = (res, j, wj)
        support.append(best[1]); w = best[2]
        for _ in range(int(rounds)):
            X = (w @ A[support]).reshape(n, 3)
            natp = O.superpose_one(nat, X).ravel()
            w = O.convex_nnls(A[support].T, natp)
        X = (w @ A[support]).reshape(n, 3)
        out[step] = (float(I.ca_rmsd(X, nat)), [int(idx[q]) for q in support], w.copy())
    return out


def oracle_weights_on(Wf, support, nat, n, natp0):
    """ORACLE. Continuous convex weights on a GIVEN support (lane O's alternating solver)."""
    r, X, w = O.oracle_hull(Wf, nat, n, support=np.asarray(support, int), natp0=natp0)
    return float(r), np.asarray(w, float)


def uniform_on(Wf, support, nat, n):
    """NATIVE-FREE construction; the RMSD readout is ORACLE (as every ceiling here is)."""
    sup = np.asarray(support, int)
    X = Wf[sup].mean(0).reshape(n, 3)
    return float(I.ca_rmsd(X, nat)), X


def fps_support(P, idx, s):
    """NATIVE-FREE farthest-point support: the medoid of `idx` under pairwise RMSD, then greedily
    the member maximising the minimum distance to those already chosen."""
    idx = np.asarray(idx, int)
    sub = P[np.ix_(idx, idx)]
    cur = [int(np.argmin(sub.mean(1)))]
    while len(cur) < min(int(s), len(idx)):
        d = sub[cur].min(0)
        d[cur] = -np.inf
        cur.append(int(np.argmax(d)))
    return [int(idx[q]) for q in cur]


# ------------------------------------------------------------------ per target
def row_for(pdb, o_row, verbose=True):
    t0 = time.time()
    cand, dis, top, order, dg = O.load_pool(pdb)
    W = cand.W; n = cand.n; nat = cand.nat_ca; rr = cand.oracle_rr        # ORACLE labels
    P = I.pairwise_rmsd(W)
    from s27 import s28_A_amp as A
    frame = A.Frame(W, np.sort(top))
    Wf = frame.Wf                                                         # (k, 3n) common frame
    natp0 = O.superpose_one(nat, frame.ref)
    natp0f = natp0.ravel()
    top75 = np.sort(top); pool = np.arange(cand.k)
    rng = np.random.default_rng(abs(hash((SALT, pdb))) % (2 ** 32))

    row: Dict = dict(pdb=pdb, n=n, fold=int(cand.fold), fail18=bool(pdb in I.FAIL18),
                     k=int(cand.k), arms={})

    def put(name, **kw):
        row["arms"][name] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                             for k, v in kw.items()}

    # ---------- anchors -------------------------------------------------------
    C, _ = O.avg_of(W, top75, P)
    put("prod", rmsd=float(I.ca_rmsd(C, nat)))                            # ORACLE readout
    # F4 reference: the argmin over the top-2**B of the score order, costing exactly B bits
    ref = {}
    for B in B_BITS:
        N = min(2 ** B, cand.k)
        ref[str(B)] = float(rr[order[:N]].min())                          # ORACLE
    row["argmin_ref"] = ref

    # ---------- the factorial: support choice x weight channel ----------------
    for Sname, idx in (("top75", top75), ("pool", pool)):
        o_items = o_row["items"]
        for s in S_FULL:
            if s > len(idx):
                continue
            key = "sparse_%s_s%d" % (Sname, s)
            sup_o = [int(q) for q in o_items[key]["support"]]             # ORACLE support, on disk
            w_o = np.asarray(o_items[key]["w"], float)
            put("A_%s_s%d" % (Sname, s), rmsd=float(o_items[key]["rmsd_cloud"]),
                sup_bits=support_bits(len(idx), s), w_bits=float("inf"))
            # B -- ORACLE support, UNIFORM weights
            rb, _ = uniform_on(Wf, sup_o, nat, n)
            put("B_%s_s%d" % (Sname, s), rmsd=rb,
                sup_bits=support_bits(len(idx), s), w_bits=0.0)
            # Q(L) -- ORACLE support, GRID weights
            for L in L_GRID:
                rq, wq = oracle_grid_weights(Wf[np.asarray(sup_o, int)], nat, n, L, w0=w_o)
                put("Q_%s_s%d_L%d" % (Sname, s, L), rmsd=rq,
                    sup_bits=support_bits(len(idx), s), w_bits=weight_bits(s, L))
            # C / D -- SCORE-PREFIX support (0 support bits)
            sup_p = [int(q) for q in order[:s]] if Sname == "pool" else [int(q) for q in order[:s]]
            rc, _ = oracle_weights_on(Wf, sup_p, nat, n, natp0f)
            rd, _ = uniform_on(Wf, sup_p, nat, n)
            put("C_%s_s%d" % (Sname, s), rmsd=rc, sup_bits=0.0, w_bits=float("inf"))
            put("D_%s_s%d" % (Sname, s), rmsd=rd, sup_bits=0.0, w_bits=0.0)
            # G / H -- DIVERSITY (farthest-point) support, NATIVE-FREE (0 support bits)
            sup_f = fps_support(P, idx, s)
            rg, _ = oracle_weights_on(Wf, sup_f, nat, n, natp0f)
            rh, _ = uniform_on(Wf, sup_f, nat, n)
            put("G_%s_s%d" % (Sname, s), rmsd=rg, sup_bits=0.0, w_bits=float("inf"))
            put("H_%s_s%d" % (Sname, s), rmsd=rh, sup_bits=0.0, w_bits=0.0)
            # E / F -- RANDOM support, MEAN over draws (min kept only to price the order statistic)
            re_, rf_ = [], []
            for _d in range(N_RAND):
                sup_r = [int(q) for q in rng.choice(idx, size=s, replace=False)]
                re_.append(oracle_weights_on(Wf, sup_r, nat, n, natp0f)[0])
                rf_.append(uniform_on(Wf, sup_r, nat, n)[0])
            put("E_%s_s%d" % (Sname, s), rmsd=float(np.mean(re_)), rmsd_min=float(np.min(re_)),
                draws=[float(x) for x in re_],
                sup_bits=support_bits(len(idx), s), w_bits=float("inf"))
            put("F_%s_s%d" % (Sname, s), rmsd=float(np.mean(rf_)), rmsd_min=float(np.min(rf_)),
                draws=[float(x) for x in rf_],
                sup_bits=support_bits(len(idx), s), w_bits=0.0)

    # ---------- the low-bit ladder: support drawn from the top-T prefix -------
    for T in T_LADDER:
        if T > cand.k:
            continue
        idxT = np.asarray(order[:T], int)
        s_max = max(s for s in S_LADDER if s <= T)
        g = oracle_greedy_support(Wf, idxT, nat, n, s_max, natp0f)        # ORACLE
        for s in S_LADDER:
            if s not in g:
                continue
            r_o, sup, w_c = g[s]
            put("T%d_s%d_cont" % (T, s), rmsd=r_o, sup_bits=support_bits(T, s),
                w_bits=float("inf"), support=[int(x) for x in sup])
            ru, _ = uniform_on(Wf, sup, nat, n)
            put("T%d_s%d_unif" % (T, s), rmsd=ru, sup_bits=support_bits(T, s), w_bits=0.0)
            for L in L_GRID:
                rq, _ = oracle_grid_weights(Wf[np.asarray(sup, int)], nat, n, L, w0=w_c)
                put("T%d_s%d_L%d" % (T, s, L), rmsd=rq, sup_bits=support_bits(T, s),
                    w_bits=weight_bits(s, L))

    # ---------- greedy-gap control: exhaustive s = 2 over the top-EXHAUST_T ---
    if cand.k >= EXHAUST_T:
        idxE = np.asarray(order[:EXHAUST_T], int)
        Ae = Wf[idxE]
        pairs = list(itertools.combinations(range(len(idxE)), 2))
        best = (np.inf, None)
        for (a, b) in pairs:                                              # ORACLE
            wj = O.convex_nnls(Ae[[a, b]].T, natp0f)
            X = (wj @ Ae[[a, b]]).reshape(n, 3)
            r = float(I.ca_rmsd(X, nat))
            if r < best[0]:
                best = (r, (a, b))
        put("EX%d_s2_cont" % EXHAUST_T, rmsd=best[0], sup_bits=support_bits(EXHAUST_T, 2),
            w_bits=float("inf"), npairs=len(pairs))

    row["secs"] = time.time() - t0
    if verbose:
        a = row["arms"]
        print("  %-6s n=%2d prod %.3f | A_pool_s2 %.3f  B(unif) %.3f  Q(L=4) %.3f | "
              "C(score) %.3f  G(fps) %.3f  E(rand) %.3f | argmin B=7 %.3f B=9 %.3f | %.1fs"
              % (pdb, n, a["prod"]["rmsd"], a["A_pool_s2"]["rmsd"], a["B_pool_s2"]["rmsd"],
                 a["Q_pool_s2_L4"]["rmsd"], a["C_pool_s2"]["rmsd"], a["G_pool_s2"]["rmsd"],
                 a["E_pool_s2"]["rmsd"], ref["7"], ref["9"], row["secs"]), flush=True)
    return row


# ------------------------------------------------------------------ phases
def jsonl_rows(f):
    out = {}
    if os.path.exists(f):
        for L in open(f):
            L = L.strip()
            if L:
                r = json.loads(L)
                out[r["pdb"]] = r
    return out


def append_row(f, row):
    os.makedirs(os.path.dirname(f), exist_ok=True)
    tmp = f + ".tmp%d" % os.getpid()
    with open(tmp, "w") as fh:
        fh.write(json.dumps(row) + "\n")
    with open(f, "a") as fh, open(tmp) as r:
        fh.write(r.read())
    os.remove(tmp)


def phase_cloud(pdbs, limit=None):
    o_rows = jsonl_rows(O_CLOUD)
    done = jsonl_rows(ROWS)
    todo = [p for p in pdbs if p not in done]
    if limit:
        todo = todo[:int(limit)]
    print("lane Q sparse ceiling: %d targets to do (%d already on disk)" % (len(todo), len(done)),
          flush=True)
    for p in todo:
        if p not in o_rows:
            print("  SKIP %s -- no lane O cloud row" % p, flush=True)
            continue
        append_row(ROWS, row_for(p, o_rows[p]))
    return jsonl_rows(ROWS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["cloud", "analyse"])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--pdbs", type=str, default=None)
    a = ap.parse_args()
    pdbs = a.pdbs.split(",") if a.pdbs else O.all_pdbs()
    if a.phase == "cloud":
        phase_cloud(pdbs, a.limit)
    else:
        from s30 import s30_Q_analyse as QA
        QA.main()


if __name__ == "__main__":
    main()
