"""SPRINT 14, OBJ -- the learning curve, against the leaked-label positive control.

Sprint 12's decisive negative was a FLAT learning curve while the identical harness with
a leaked label was already loud at n=8.  This runs that diagnostic on the torsion-space
objective.

It is cheap because normal equations are ADDITIVE over per-target centered blocks: one
featurization pass over the training pool gives a per-target block, and training on ANY
subset of m targets is just the sum of m blocks.  No refitting, no re-featurizing.
Blocks are stored float32 (68 MB each) to fit the box.

Test set is one whole pinned fold; the pool is every target in the OTHER folds, so no
training subset ever contains a fold-mate of a test target.

    python -m s14.obj_lc
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                    # noqa: E402
from s13.qarch_lib import spearman                 # noqa: E402
from s14 import obj_enum as E                      # noqa: E402
from s14 import obj_model as M                     # noqa: E402
from s14 import obj_train as T                     # noqa: E402
from s14 import obj_exp as X                       # noqa: E402
from s14 import obj_headline as H                  # noqa: E402

SIZES = [1, 2, 4, 8]
NREP = 4


def target_block(fz, p, nuni, seed):
    """One target's centered normal-equation block, stored float32."""
    en = T.Enum(p)
    idx = T.sample_idx(en, nuni, 0, np.random.default_rng(seed + hash(p) % 9973))
    Xm = fz.csr(p, T._decode(idx, en.n, en.k))
    y = en.pct[idx].astype(np.float64)
    N = Xm.shape[0]
    mu = np.asarray(Xm.sum(0)).ravel() / N
    ybar = float(y.mean())
    XtX = (Xm.T @ Xm).toarray() - N * np.outer(mu, mu)
    Xty = Xm.T @ y - N * mu * ybar
    yty = float(y @ y) - N * ybar * ybar
    return XtX.astype(np.float32), Xty.astype(np.float64), yty, N


def solve_sum(blocks, lam, leak=False):
    XtX = np.zeros((M.DIM, M.DIM)); Xty = np.zeros(M.DIM); yty = 0.0; n = 0
    for a, b, c, N in blocks:
        XtX += a; Xty += b; yty += c; n += N
    if leak:
        d = M.DIM
        A = np.zeros((d + 1, d + 1)); bb = np.zeros(d + 1)
        A[:d, :d] = XtX; A[:d, d] = A[d, :d] = Xty; A[d, d] = yty
        bb[:d] = Xty; bb[d] = yty
        XtX, Xty = A, bb
    A = XtX / max(n, 1); b = Xty / max(n, 1)
    s = np.sqrt(np.maximum(np.diag(A), 0.0)); s[s < 1e-9] = 1.0
    D = 1.0 / s
    A = A * D[:, None] * D[None, :] + H.LAM * np.eye(len(A))
    return np.linalg.solve(A, b * D) * D


def run(targets, seed=0, nuni=8000):
    fz = H.fz_for(True)
    folds = {p: T.Enum(p).fold for p in targets}
    tf = max(set(folds.values()), key=lambda f: sum(1 for p in targets if folds[p] == f))
    test = [p for p in targets if folds[p] == tf]
    pool = [p for p in targets if folds[p] != tf]
    print(f"held-out fold {tf} ({len(test)}): {test}")
    print(f"pool ({len(pool)}): {pool}", flush=True)
    evx = {p: X.eval_idx(T.Enum(p)) for p in test}
    rr = {p: T.Enum(p).rmsd[evx[p]].astype(np.float64) for p in test}
    Sdec = {p: T._decode(evx[p], T.Enum(p).n, T.Enum(p).k) for p in test}

    t0 = time.time()
    blk = {}
    for p in pool:
        blk[p] = target_block(fz, p, nuni, seed)
        print(f"  block {p} [{time.time()-t0:.0f}s]", flush=True)

    sizes = [m for m in SIZES if m <= len(pool)]
    if len(pool) not in sizes:
        sizes.append(len(pool))
    rows = {}
    for arm, leak in (("learned", False), ("leak_control", True)):
        for m in sizes:
            acc = []
            reps = NREP if m < len(pool) else 1
            for rep in range(reps):
                rng = np.random.default_rng(7919 * m + rep + seed)
                tr = list(rng.choice(pool, m, replace=False))
                w = solve_sum([blk[q] for q in tr], H.LAM, leak=leak)
                for p in test:
                    e = fz.score_with(np.asarray(w)[:M.DIM], p, Sdec[p])
                    if leak:
                        e = e + float(w[M.DIM]) * T.Enum(p).pct[evx[p]].astype(np.float64)
                    r = rr[p]
                    o = np.argsort(e, kind="mergesort")
                    dec = o[:max(50, len(e) // 10)]
                    acc.append(dict(d_decile=float(r[dec].mean() - r.mean()),
                                    d_top100=float(r[o[:100]].mean() - r.mean()),
                                    rho_global=spearman(e, r),
                                    rho_decile=spearman(e[dec], r[dec])))
            v = {k: float(np.nanmean([a[k] for a in acc])) for k in acc[0]}
            v["n_obs"] = len(acc)
            v["sd_d_decile"] = float(np.std([a["d_decile"] for a in acc]))
            rows[f"{arm}|{m}"] = v
            print(f"{arm:13s} m={m:2d}  d_decile={v['d_decile']:+.3f} "
                  f"(sd {v['sd_d_decile']:.3f})  d_top100={v['d_top100']:+.3f}  "
                  f"rho_g={v['rho_global']:+.3f}  rho_dec={v['rho_decile']:+.3f}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    I.write("s14_obj_lc", {"test_fold": int(tf), "test": test, "pool": pool,
                           "n_uniform": nuni, "lam": H.LAM, "rows": rows},
            n_expected=1)
    return rows


if __name__ == "__main__":
    tg = [a for a in sys.argv[1:] if not a.startswith("-")]
    run(tg or [t["pdb"] for t in I.targets() if E.have(t["pdb"])])
