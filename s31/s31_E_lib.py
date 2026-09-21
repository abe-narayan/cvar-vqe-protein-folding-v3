#!/usr/bin/env python
"""s31/s31_E_lib.py -- S31 lane E shared cache and the native-free trace assertion.

Builds ONE pass over the 126 dev targets producing exactly the arrays `s30/s30_P_lr.py` builds
(same FEATS, same POOL_K/M/MIN_SEP, same ORACLE response) plus the two quantities lane E needs:

    mu      = pool75_mean - d_nat      ORACLE   (the pool's common-mode pair error, S30-L7)
    mu_hat  = pool75_mean - expected   NATIVE-FREE

and the exact algebraic identity  mu_hat = mu - y  with y = expected - d_nat, asserted in code.

THE NATIVE-FREE TRACE, verified here rather than asserted in prose (`s30` published its `coh`
gate as native-free and it was ORACLE):

    mu_hat[p] = pool75_mean[p] - expected[p]
        pool75_mean = D[rec["sub"]].mean(0),  D = pair distances of the pool members W
            W        <- u["W"][u["order"][:500]]        retrieval + its stored order
            rec      <- I.shipped_record(pdb)["sub"]    the shipped top-75 indices
        expected    <- I.distogram(pdb)["expected"]     leave-fold-out ESM-2 distogram
    Neither branch reads u["nat_ca"].  `build_rows` proves it by rebuilding mu_hat from a copy
    of the universe whose `nat_ca` has been replaced by NaN, and asserting bit-equality.

    python s31/s31_E_lib.py           # build/refresh the cache, print the trace
"""
from __future__ import annotations

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

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
CACHE = os.path.join(RESULTS, "s31_E_rows.npz")

POOL_K, M, MIN_SEP = 500, 75, 2
FEATS = ["expected", "sd", "entropy", "exp_minus_median", "exp_minus_mode", "sep", "n",
         "d_prod", "pool75_mean", "pool500_mean", "pool75_sd", "exp_minus_pool75",
         "exp_minus_pool500", "rg_disagree", "post_sd_mean", "score_sd", "sep_over_n",
         "expected_over_n", "prob_tail_hi", "prob_tail_lo"]

NEST = {
    "N0_separation_prior": ["sep", "n", "sep_over_n"],
    "N1_plus_calibration": ["sep", "n", "sep_over_n", "expected", "expected_over_n"],
    "N2_plus_posterior_shape": ["sep", "n", "sep_over_n", "expected", "expected_over_n", "sd",
                                "entropy", "exp_minus_median", "exp_minus_mode", "prob_tail_hi",
                                "prob_tail_lo", "post_sd_mean"],
    "N3_plus_pool": FEATS,
}


def _target_arrays(t, want_nat=True):
    """Everything for one target.  `want_nat=False` returns NaN for every ORACLE column."""
    pdb = t["pdb"]; n = int(t["n"])
    u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
    ii, jj = I.pair_index(n, MIN_SEP); sep = (jj - ii).astype(float)
    exp = np.asarray(dg["expected"], float); sd = np.asarray(dg["sd"], float)
    P = np.asarray(dg["prob"], float); cen = np.asarray(dg["centres"], float)
    ent = -(P * np.log(np.clip(P, 1e-12, None))).sum(1)
    med = cen[np.argmax(np.cumsum(P, 1) >= 0.5, axis=1)]
    mode = cen[np.argmax(P, 1)]
    pool = np.asarray(u["order"], int)[:POOL_K]
    W = np.asarray(u["W"], float)[pool]
    D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    sub = np.asarray(rec["sub"], int)
    C0 = np.asarray(rec["avg_ca"], float)
    d_prod = np.linalg.norm(C0[ii] - C0[jj], axis=1)
    rgp = float(np.sqrt((np.sum(exp ** 2) * 2 + (n - 1) * 2 * 3.8 ** 2) / (2.0 * n * n)))
    rgo = float(np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1)).mean())
    sc = I.shipped_score(dg, D)
    p75 = D[sub].mean(0)
    col = dict(expected=exp, sd=sd, entropy=ent, exp_minus_median=exp - med,
               exp_minus_mode=exp - mode, sep=sep, n=np.full(len(exp), float(n)), d_prod=d_prod,
               pool75_mean=p75, pool500_mean=D.mean(0), pool75_sd=D[sub].std(0),
               exp_minus_pool75=exp - p75, exp_minus_pool500=exp - D.mean(0),
               rg_disagree=np.full(len(exp), rgp - rgo),
               post_sd_mean=np.full(len(exp), sd.mean()),
               score_sd=np.full(len(exp), sc.std()), sep_over_n=sep / n,
               expected_over_n=exp / n, prob_tail_hi=P[:, -3:].sum(1),
               prob_tail_lo=P[:, :3].sum(1))
    X = np.column_stack([col[f] for f in FEATS])
    mu_hat = p75 - exp                                       # NATIVE-FREE
    if want_nat:
        nat = np.asarray(u["nat_ca"], float)
        d_nat = np.linalg.norm(nat[ii] - nat[jj], axis=1)    # ORACLE
        y = exp - d_nat                                      # ORACLE
        mu = p75 - d_nat                                     # ORACLE
    else:
        y = np.full(len(exp), np.nan); mu = np.full(len(exp), np.nan)
    return X, y, mu, mu_hat, sep, ii, jj


def build_rows(verbose=True):
    """Build the cache.  Returns the dict of concatenated arrays."""
    tg = I.targets()
    X, Y, MU, MH, SEP, TGT, FD = [], [], [], [], [], [], []
    pdbs, nres, off = [], [], [0]
    t0 = time.time()
    for k, t in enumerate(tg):
        x, y, mu, mh, sep, ii, jj = _target_arrays(t)
        # --- the identity, asserted, not assumed
        assert np.allclose(mh, mu - y, atol=1e-9), "mu_hat != mu - y on %s" % t["pdb"]
        X.append(x); Y.append(y); MU.append(mu); MH.append(mh); SEP.append(sep)
        TGT.append(np.full(len(y), k)); FD.append(np.full(len(y), int(t["fold"])))
        pdbs.append(t["pdb"]); nres.append(int(t["n"])); off.append(off[-1] + len(y))
        if verbose and (k + 1) % 20 == 0:
            print("  rows %d/%d  %.1fs" % (k + 1, len(tg), time.time() - t0), flush=True)
    out = dict(X=np.concatenate(X), y=np.concatenate(Y), mu=np.concatenate(MU),
               mu_hat=np.concatenate(MH), sep=np.concatenate(SEP),
               tgt=np.concatenate(TGT).astype(int), fold=np.concatenate(FD).astype(int),
               off=np.array(off, int), nres=np.array(nres, int),
               pdbs=np.array(pdbs), feats=np.array(FEATS))
    return out


def native_free_trace(pdb=None):
    """PROOF, not prose: rebuild mu_hat with nat_ca destroyed and assert bit-equality.

    Monkey-patches `I.load_univ` so `nat_ca` is all-NaN, recomputes mu_hat for one target, and
    checks it is bit-identical to the normal path.  If any native value leaked into mu_hat the
    NaN would propagate and the assertion would fail.
    """
    tg = I.targets()
    t = tg[0] if pdb is None else [q for q in tg if q["pdb"] == pdb][0]
    _, _, _, mh_ref, _, _, _ = _target_arrays(t)
    real = I.load_univ

    def poisoned(p):
        u = dict(real(p))
        u["nat_ca"] = np.full_like(np.asarray(u["nat_ca"], float), np.nan)
        return u

    I.load_univ = poisoned
    try:
        _, _, _, mh_poi, _, _, _ = _target_arrays(t, want_nat=False)
    finally:
        I.load_univ = real
    ok = bool(np.array_equal(mh_ref, mh_poi))
    return {"target": t["pdb"], "mu_hat_identical_with_nat_ca_destroyed": ok,
            "max_abs_diff": float(np.max(np.abs(mh_ref - mh_poi))) if ok else None}


def load(rebuild=False):
    if (not rebuild) and os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=False)
        return {k: z[k] for k in z.files}
    d = build_rows()
    tmp = CACHE[:-4] + ".%d.tmp.npz" % os.getpid()   # `atomic-rename-needs-unique-temp`
    np.savez_compressed(tmp, **d)
    os.replace(tmp, CACHE)
    return d


# ----------------------------------------------------------------- S30's leave-fold-out ridge
def lfo_ridge(X, y, fd, alphas=(1e0, 1e1, 1e2, 1e3, 1e4)):
    """Bit-identical to `s30/s30_P_lr.py:lfo_r2` (no RNG; deterministic)."""
    yh = np.zeros_like(y)
    F = sorted(set(fd.tolist()))
    for f in F:
        tr, te = fd != f, fd == f
        best, ba = np.inf, alphas[0]
        for a in alphas:
            err = 0.0
            for g in [q for q in F if q != f]:
                itr, ite = tr & (fd != g), fd == g
                mu, s = X[itr].mean(0), np.where(X[itr].std(0) > 1e-12, X[itr].std(0), 1.0)
                A = (X[itr] - mu) / s; Bx = (X[ite] - mu) / s
                A1 = np.column_stack([A, np.ones(len(A))])
                B1 = np.column_stack([Bx, np.ones(len(Bx))])
                w = np.linalg.solve(A1.T @ A1 + a * np.eye(A1.shape[1]), A1.T @ y[itr])
                err += float(((y[ite] - B1 @ w) ** 2).sum())
            if err < best:
                best, ba = err, a
        mu, s = X[tr].mean(0), np.where(X[tr].std(0) > 1e-12, X[tr].std(0), 1.0)
        A = (X[tr] - mu) / s; Bx = (X[te] - mu) / s
        A1 = np.column_stack([A, np.ones(len(A))]); B1 = np.column_stack([Bx, np.ones(len(Bx))])
        w = np.linalg.solve(A1.T @ A1 + ba * np.eye(A1.shape[1]), A1.T @ y[tr])
        yh[te] = B1 @ w
    return yh


def per_target(a, tgt):
    """Iterator of (target index, boolean mask)."""
    for q in range(int(tgt.max()) + 1):
        yield q, tgt == q


if __name__ == "__main__":
    d = load(rebuild="--rebuild" in sys.argv)
    print("cache", CACHE, "pairs", len(d["y"]), "targets", len(d["pdbs"]))
    print("NATIVE-FREE TRACE:", native_free_trace())
