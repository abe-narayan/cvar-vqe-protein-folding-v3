"""SPRINT 18 / MATH -- shared machinery for the functional-ANOVA derivation.

Nothing here is an experiment.  It is the algebra the two experiment modules
(`s18/math_lattice.py`, `s18/math_anova.py`) share, plus the artefact conventions the brief
mandates (deterministic config hash, completion flag, persisted coefficients).

THE ONE DEFINITION EVERYTHING RESTS ON.  Let `X = prod_i X_i` be a product space and
`mu = tensor_i mu_i` a PRODUCT probability measure on it.  For `E in L^2(mu)` the
first-order functional ANOVA decomposition is

    E_0     = E_mu[E]
    f_i(t)  = E_mu[E | theta_i = t] - E_0                       (so E_{mu_i}[f_i] = 0)
    E_le1   = E_0 + sum_i f_i(theta_i)
    E_ge2   = E - E_le1

`E_le1` is the ORTHOGONAL PROJECTION of `E` onto the closed subspace
`V_1 = {c + sum_i g_i(theta_i)}` of mu-square-integrable additive functions.  Orthogonality
is what makes the split unique: `<E_ge2, g> = 0` for every `g in V_1`, and hence
`Var_mu(E) = Var_mu(E_le1) + Var_mu(E_ge2)`.  Both are checked numerically in
`s18/math_lattice.py`; neither is assumed.

TWO FACTORISATIONS, TWO DIFFERENT OBJECTS.  The same energy admits an ANOVA under any
product structure.  On the enumerated lattice there are two natural ones:

  * QUBIT factorisation, `X = {0,1}^{n_q}`.  Order-<=1 ANOVA under the uniform measure is
    EXACTLY the Walsh weight-<=1 projection (`P_le1_qubit`).
  * RESIDUE factorisation, `X = {0..k-1}^n`.  Order-<=1 ANOVA is the residue-additive field
    (`P_le1_res`), which in Walsh language is weight-0 + weight-1 + ALL INTRA-RESIDUE
    weight-2 coefficients.

`P_le1_qubit` is a STRICT subspace of `P_le1_res` when k > 2.  The gap is not a detail: it
is exactly the part of the per-residue field that a binary encoding hides above degree 1.

Everything here is native-free.  No function in this file reads an RMSD, a native structure
or any label.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(HERE, "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

SALT = "s18math"


# ===================================================================== artefacts
def cfg_hash(d):
    """Deterministic hash of a config dict -- persisted beside every artefact."""
    s = json.dumps(d, sort_keys=True, default=str)
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def _enc(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


_CK = {}


def ck(tag, key, value):
    """Merge-on-write incremental checkpoint into `s18/results/math_<tag>.json`.

    `_complete` is written ONLY by `seal`.  The brief's rule -- never treat a partial file as
    complete -- is enforced by every reader checking that flag.
    """
    path = os.path.join(RESULTS, f"math_{tag}.json")
    if os.path.exists(path):
        try:
            with open(path) as fh:
                disk = json.load(fh)
            mem = _CK.setdefault(tag, {})
            for k, v in disk.items():
                mem.setdefault(k, v)
        except Exception:
            pass
    d = _CK.setdefault(tag, {})
    d[key] = value
    d["_written"] = time.strftime("%Y-%m-%d %H:%M:%S")
    d.setdefault("_complete", False)
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(d, fh, indent=1, default=_enc)
    for _ in range(30):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.4)
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1, default=_enc)


def ck_load(tag):
    path = os.path.join(RESULTS, f"math_{tag}.json")
    if os.path.exists(path):
        with open(path) as fh:
            d = json.load(fh)
        _CK[tag] = d
        return d
    return _CK.setdefault(tag, {})


def seal(tag, config):
    d = ck_load(tag)
    d["_config"] = config
    d["_config_hash"] = cfg_hash(config)
    d["_complete"] = True
    d["_sealed"] = time.strftime("%Y-%m-%d %H:%M:%S")
    path = os.path.join(RESULTS, f"math_{tag}.json")
    tmp = path + f".seal{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(d, fh, indent=1, default=_enc)
    os.replace(tmp, path)
    return d["_config_hash"]


# ===================================================== Walsh (uniform Boolean cube)
def fwht(a):
    """Normalised fast Walsh-Hadamard transform: `c = FWHT(f)/2^n` are the coefficients of
    `f = sum_S c_S chi_S`, `chi_S(x) = prod_{q in S} (-1)^{x_q}`.  Byte-identical in
    convention to `s17/q_theory.py` so the reproduction is exact."""
    a = np.array(a, dtype=np.float64, copy=True)
    n = a.size
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(-1)
        h *= 2
    return a / n


def ifwht(c):
    a = np.array(c, dtype=np.float64, copy=True)
    n = a.size
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(-1)
        h *= 2
    return a


def popcount(a):
    a = np.asarray(a, np.int64)
    c = np.zeros_like(a)
    x = a.copy()
    while x.any():
        c += (x & 1)
        x >>= 1
    return c


def walsh_truncate(f, nq, d):
    """Zero every Walsh coefficient of weight > d and invert.  `d = 1` is the object
    Sprint 17 called `degree-1`."""
    c = fwht(f)
    deg = popcount(np.arange(1 << nq, dtype=np.int64))
    c = c.copy()
    c[deg > d] = 0.0
    return ifwht(c)


def walsh_residue_additive(f, n, bits_per_res):
    """Keep weight-0, weight-1 and every weight-2 coefficient whose two qubits lie in the
    SAME residue block; zero everything else.  Claimed EQUAL to the residue-ANOVA order-1
    projection under uniform mu; `math_lattice.py` verifies that to machine precision."""
    nq = n * bits_per_res
    idx = np.arange(1 << nq, dtype=np.int64)
    deg = popcount(idx)
    # qubit q (0 = MSB) belongs to residue q // bits_per_res
    keep = deg <= 1
    if bits_per_res > 1:
        two = np.flatnonzero(deg == 2)
        b = idx[two]
        hi = nq - 1 - np.floor(np.log2(b)).astype(np.int64)
        lo = nq - 1 - np.floor(np.log2(b - (np.int64(1) << (nq - 1 - hi)))).astype(np.int64)
        same = (hi // bits_per_res) == (lo // bits_per_res)
        keep[two[same]] = True
    c = fwht(f)
    c = np.where(keep, c, 0.0)
    return ifwht(c)


# ============================================ exact ANOVA on a discrete product space
def cond_means(T, p, i):
    """`E_mu[E | x_i = s]` for every `s`, EXACT.

    `T` is the energy as a `(k,)*n` tensor (axis a = residue a), `p` is `(n, k)` with rows
    summing to 1.  Contracting the axes from the last downward keeps the index of every
    axis that has not yet been contracted, so no bookkeeping is needed.
    """
    n = T.ndim
    A = T
    for a in reversed(range(n)):
        if a == i:
            continue
        A = np.tensordot(A, p[a], axes=([a], [0]))
    return np.asarray(A, float)


def anova1_discrete(E, n, k, p=None):
    """Order-<=1 functional ANOVA of a fully enumerated `E` under the product measure `p`.

    Returns `(E0, f, E_le1)` with `f` shaped `(n, k)` and `sum_s p[i,s] f[i,s] = 0` for
    every `i` -- the centring that makes the decomposition orthogonal.
    """
    E = np.asarray(E, float)
    if p is None:
        p = np.full((n, k), 1.0 / k)
    p = np.asarray(p, float)
    assert p.shape == (n, k)
    assert np.allclose(p.sum(1), 1.0)
    T = E.reshape((k,) * n)
    # E0 by full contraction against the product measure
    A = T
    for a in reversed(range(n)):
        A = np.tensordot(A, p[a], axes=([a], [0]))
    E0 = float(A)
    f = np.empty((n, k), float)
    for i in range(n):
        m = cond_means(T, p, i)
        f[i] = m - E0
        f[i] -= float(p[i] @ f[i])          # exact centring; removes float drift
    # reconstruct on the full lattice
    pw = k ** np.arange(n - 1, -1, -1)
    idx = np.arange(k ** n, dtype=np.int64)
    S = (idx[:, None] // pw[None, :]) % k
    E_le1 = E0 + f[np.arange(n)[None, :], S].sum(1)
    return E0, f, E_le1


def states_of(N, n, k):
    pw = k ** np.arange(n - 1, -1, -1)
    idx = np.arange(N, dtype=np.int64)
    return (idx[:, None] // pw[None, :]) % k


# ===================================================================== statistics
def paired(a, b, folds=None, n_boot=4000, seed=0):
    """Paired a - b with a fold-CLUSTERED bootstrap when folds are given.  Negative = a
    better.  Median, W/L and the drop-top diagnostics beside the mean, per standing law."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = a - b
    n = len(d)
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    ci_cl = None
    if folds is not None:
        folds = np.asarray(folds)
        groups = [np.flatnonzero(folds == f) for f in np.unique(folds)]
        bc = np.empty(n_boot)
        for t in range(n_boot):
            pick = rng.integers(0, len(groups), len(groups))
            sel = np.concatenate([groups[q] for q in pick])
            bc[t] = d[sel].mean()
        ci_cl = [float(np.percentile(bc, 2.5)), float(np.percentile(bc, 97.5))]
    order = np.argsort(d)
    out = {"n": int(n), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
           "mean_diff": float(d.mean()), "median_diff": float(np.median(d)),
           "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
           "n_better": int((d < 0).sum()), "n_worse": int((d > 0).sum()),
           "n_tie": int((d == 0).sum()),
           "drop_top3_mean_diff": float(d[order[3:]].mean()) if n > 3 else None}
    if folds is not None:
        out["ci95_foldclustered"] = ci_cl
        pf = {int(f): float(d[folds == f].mean()) for f in np.unique(folds)}
        out["per_fold"] = pf
        sgn = np.sign(d.mean())
        out["folds_same_sign"] = f"{sum(1 for v in pf.values() if np.sign(v) == sgn)}" \
                                 f"/{len(pf)}"
    return out


def spearman(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / den) if den > 0 else 0.0
