#!/usr/bin/env python
"""s30/s30_Q_quadric.py -- lane Q: what the QUADRIC (second-moment) reachable-tail class can reach,
measured BEFORE any Hamiltonian is built.

The coordinator's construction: with an endogenous order the cost is c_x = <grad f(R_lambda), W_x>,
so candidates are sorted by projection onto ONE self-consistent direction and the reachable tails
are HALFSPACE cuts. Adding a second-moment term V = f(sum l W, sum l W W^T) makes grad V_x
QUADRATIC in W_x, so the cuts become QUADRICS.

This measures the ORACLE ceiling of each reachable-tail class on the real pools:

    EXOGENOUS   prefixes of the deployed DIS order                  (what ships)
    LINEAR      prefixes of the order induced by <c, W_x>           (halfspace cuts)
    QUADRIC     prefixes of the order induced by W_x^T A W_x + <c, W_x>
    FREE        the unreachable per-size optimum, for scale

Each class is read out the SAME way -- the uniform coordinate average of the selected prefix --
so the comparison isolates WHICH SET, which is the only thing the class changes.

ORDER-STATISTIC DISCIPLINE (contract rule 8): the ceiling over K sampled directions is a
best-of-K. The ceiling is therefore reported as a CURVE in K, priced with `best_of_k_within`, and
separated from the NATIVE-FREE RULE arms (PC1, medoid distance, Rg, the DIS score), which are
single directions produced by a rule and are the only deployable rows.

EVERY ceiling here is ORACLE. No deployable parameter is chosen.

    python s30/s30_Q_quadric.py [--limit N] [--ndir K]
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s29 import s29_O_ladder as O          # noqa: E402

ROWS = os.path.join(HERE, "results", "s30_Q_quadric_rows.jsonl")
OUT = os.path.join(HERE, "results", "s30_Q_quadric.json")
SALT = "s30Qquad"
NDIR = 128                 # sampled directions per class
K_CURVE = (1, 2, 4, 8, 16, 32, 64, 128)
CHUNK = 16                 # directions per Kabsch batch


def prefix_rmsds(Wf, orderings, nat, n):
    """For each ordering (D, k) return the RMSD of the uniform average of every prefix m = 1..k.
    Returns (D, k). The frame is the deployed shared frame, so a prefix average is a running mean."""
    D, k = orderings.shape
    out = np.empty((D, k))
    for a in range(0, D, CHUNK):
        b = min(a + CHUNK, D)
        idx = orderings[a:b]                                  # (d, k)
        Z = Wf[idx]                                           # (d, k, 3n)
        run = np.cumsum(Z, axis=1) / np.arange(1, k + 1)[None, :, None]
        X = run.reshape(-1, n, 3)
        out[a:b] = I.kabsch_rmsd_batch(X, nat).reshape(b - a, k)   # ORACLE
    return out


def row_for(pdb, ndir=NDIR, verbose=True):
    t0 = time.time()
    cand, dis, top, order, dg = O.load_pool(pdb)
    W = cand.W; n = cand.n; nat = cand.nat_ca; k = int(cand.k)     # ORACLE label: nat
    from s27 import s28_A_amp as A
    frame = A.Frame(W, np.sort(top))
    Wf = frame.Wf                                                   # (k, 3n), deployed frame
    Wc = Wf - Wf.mean(0, keepdims=True)
    d = Wf.shape[1]
    rng = np.random.default_rng(abs(hash((SALT, pdb))) % (2 ** 32))
    row = dict(pdb=pdb, n=n, k=k, d=int(d), fold=int(cand.fold),
               fail18=bool(pdb in I.FAIL18), classes={})

    # ---- EXOGENOUS: the deployed DIS order (1 ordering)
    exo = prefix_rmsds(Wf, np.asarray(order, int)[None, :], nat, n)[0]
    row["classes"]["exogenous"] = dict(curve=exo.tolist(), best=float(exo.min()),
                                       argbest_m=int(exo.argmin()) + 1,
                                       at_m75=float(exo[74]))

    # ---- FREE: the unreachable per-size optimum, for scale.  Greedy forward selection on the
    #      running mean is NOT the exact optimum; the exact best-of-size-1 and the ORACLE sort
    #      (candidates ordered by their OWN rmsd) bound it.
    rr = cand.oracle_rr                                             # ORACLE
    oracle_sort = np.argsort(rr).astype(int)
    fre = prefix_rmsds(Wf, oracle_sort[None, :], nat, n)[0]
    row["classes"]["oracle_sort"] = dict(curve=fre.tolist(), best=float(fre.min()),
                                         argbest_m=int(fre.argmin()) + 1)

    # ---- LINEAR: halfspace cuts, ndir random directions
    C = rng.normal(size=(ndir, d))
    C /= np.linalg.norm(C, axis=1, keepdims=True)
    sl = C @ Wc.T                                                   # (ndir, k)
    ordl = np.argsort(sl, axis=1)
    rl = prefix_rmsds(Wf, ordl, nat, n)                             # (ndir, k)
    bl = rl.min(1)                                                  # best prefix per direction
    ml = rl.argmin(1) + 1                                           # WHERE the best prefix sits
    row["classes"]["linear"] = dict(per_dir_best=bl.tolist(), per_dir_m=ml.tolist(),
                                    best=float(bl.min()), best_m=int(ml[int(bl.argmin())]),
                                    at_m75=float(rl[:, 74].min()),
                                    at_m75_curveK=[float(rl[:K, 74].min()) for K in K_CURVE if K <= ndir],
                                    frac_m_le_2=float((ml <= 2).mean()),
                                    curve_bestK=[float(bl[:K].min()) for K in K_CURVE if K <= ndir])

    # ---- QUADRIC: second-moment cuts. A = sum of r rank-1 terms (r = 3), plus a linear part.
    rq = []
    for _ in range(ndir):
        U = rng.normal(size=(3, d))
        U /= np.linalg.norm(U, axis=1, keepdims=True)
        sgn = rng.choice([-1.0, 1.0], size=3)
        q = ((Wc @ U.T) ** 2 * sgn[None, :]).sum(1)                 # x^T A x, A = sum sgn_i u_i u_i^T
        c = rng.normal(size=d); c /= np.linalg.norm(c)
        rq.append(q + (Wc @ c))
    sq = np.asarray(rq)                                             # (ndir, k)
    ordq = np.argsort(sq, axis=1)
    rqq = prefix_rmsds(Wf, ordq, nat, n)
    bq = rqq.min(1)
    mq = rqq.argmin(1) + 1
    row["classes"]["quadric"] = dict(per_dir_best=bq.tolist(), per_dir_m=mq.tolist(),
                                     best=float(bq.min()), best_m=int(mq[int(bq.argmin())]),
                                     at_m75=float(rqq[:, 74].min()),
                                     at_m75_curveK=[float(rqq[:K, 74].min()) for K in K_CURVE if K <= ndir],
                                     frac_m_le_2=float((mq <= 2).mean()),
                                     curve_bestK=[float(bq[:K].min()) for K in K_CURVE if K <= ndir])

    # ---- NATIVE-FREE RULE directions: the only DEPLOYABLE rows in this file.
    #      Each is a rule that produces one direction/score from the pool without the native.
    P = I.pairwise_rmsd(W)
    med = int(I.medoid(P))
    U_, S_, Vt_ = np.linalg.svd(Wc, full_matrices=False)
    rules = {}
    for p in range(3):
        rules["PC%d" % (p + 1)] = Wc @ Vt_[p]
    rules["dist_to_medoid"] = P[med]
    rules["rg"] = np.linalg.norm(W - W.mean(1, keepdims=True), axis=2).mean(1)
    rules["dis_score"] = np.asarray(dis, float)
    # the second-moment rule the coordinator's construction points at: the pool's own dispersion,
    # i.e. the squared deviation of each candidate from the pool mean (a QUADRATIC form in W_x)
    rules["disp2"] = (Wc ** 2).sum(1)
    names = sorted(rules)
    ordr = np.asarray([np.argsort(rules[nm]) for nm in names], int)
    ordr = np.concatenate([ordr, ordr[:, ::-1]], axis=0)            # both sign conventions
    rr_rule = prefix_rmsds(Wf, ordr, nat, n)
    row["rules"] = {}
    for a, nm in enumerate(names + [x + "_rev" for x in names]):
        cur = rr_rule[a]
        row["rules"][nm] = dict(best=float(cur.min()), argbest_m=int(cur.argmin()) + 1,
                                at_m75=float(cur[74]), curve=cur.tolist())

    row["secs"] = time.time() - t0
    if verbose:
        cl, cq = row["classes"]["linear"], row["classes"]["quadric"]
        print("  %-6s n=%2d | exo %.3f (m=%d) | LIN %.3f (m=%d, %.0f%% dirs peak at m<=2; @m75 %.3f) | "
              "QUAD %.3f (m=%d; @m75 %.3f) | oracle-sort %.3f | %.1fs"
              % (pdb, n, row["classes"]["exogenous"]["best"], row["classes"]["exogenous"]["argbest_m"],
                 cl["best"], cl["best_m"], 100 * cl["frac_m_le_2"], cl["at_m75"],
                 cq["best"], cq["best_m"], cq["at_m75"],
                 row["classes"]["oracle_sort"]["best"], row["secs"]), flush=True)
    return row


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
    with open(f, "a") as fh:
        fh.write(json.dumps(row) + "\n")


def main(limit=None, ndir=NDIR):
    pdbs = O.all_pdbs()
    done = jsonl_rows(ROWS)
    todo = [p for p in pdbs if p not in done]
    if limit:
        todo = todo[:int(limit)]
    print("quadric reachable-class ceiling: %d to do (%d done), %d directions per class"
          % (len(todo), len(done), ndir), flush=True)
    for p in todo:
        append_row(ROWS, row_for(p, ndir))
    return jsonl_rows(ROWS)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--ndir", type=int, default=NDIR)
    a = ap.parse_args()
    main(a.limit, a.ndir)
