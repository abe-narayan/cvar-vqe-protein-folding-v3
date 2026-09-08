"""SPRINT 13 -- PAULI-SPECTRUM: shared machinery.

THE OBJECT.  A peptide conformation in this architecture is a bitstring: each residue picks
one of k discrete torsion states (`torsion_lib2.library_for`), encoded in b = log2(k) qubits,
so an L-residue peptide is m = L*b qubits.  Both energy models -- the Legacy 11-term
knowledge-based potential and genuine AMBER ff14SB/GBn2 -- are FUNCTIONS OF THE BITSTRING,
hence DIAGONAL in the computational basis, hence exactly a weighted sum of Pauli-Z strings

    E(z) = sum_S c_S * prod_{i in S} z_i,        z_i = +-1,  S subset of {0..m-1}

which is the Walsh-Hadamard expansion.  `c_S` is exact and computable by a fast
Walsh-Hadamard transform over the full 2^m energy table, and estimable above that size by
the unbiased Krawtchouk pair estimator below.

DEFINITIONS USED THROUGHOUT (all on the uniform product measure over configurations):

    total variance      Var = sum_{S != empty} c_S^2      (= E[f^2] - E[f]^2, Parseval)
    weight share        V_w = sum_{|S| = w} c_S^2 / Var
    mean Pauli weight   W_mean = sum_w w * V_w
    effective weight    W_eff  = 1 / sum_w V_w^2          (participation number over w)
    tail                T_{>3} = sum_{w>3} V_w

and, encoding-INDEPENDENTLY, the residue-level functional ANOVA order

    D_d = sum_{S : |res(S)| = d} c_S^2 / Var

where res(S) is the set of residues owning at least one bit of S.  D_d is a property of the
energy on the k-ary configuration space and NO index encoding can change it; the encoding
can only redistribute each D_d over Pauli weights w in [d, d*b].  That inequality is
asserted in `check_order_weight_bound`.

NOTHING here reads native coordinates.
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s13 import qarch_lib as Q            # noqa: E402  (state space, energies -- reused)

RESULTS = os.path.join(ROOT, "s13", "results")
CACHE = os.path.join(ROOT, "s13", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)


# ============================================================ the transform
def fwht(a):
    """In-place-safe fast Walsh-Hadamard transform.  Returns  A_S = sum_x (-1)^{S.x} a_x."""
    a = np.array(a, dtype=float, copy=True)
    n = a.size
    assert n & (n - 1) == 0, "length must be a power of two"
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(n)
        h *= 2
    return a


def walsh_coeffs(table):
    """Exact Pauli-Z coefficients c_S of the diagonal operator given by `table` (len 2^m)."""
    t = np.asarray(table, float).ravel()
    return fwht(t) / t.size


def ifwht_reconstruct(c):
    """Inverse: rebuild the energy table from its coefficients (for the exactness check)."""
    return fwht(np.asarray(c, float))


POPC = np.unpackbits(np.arange(1 << 16, dtype=np.uint16).view(np.uint8)
                     ).reshape(-1, 16).sum(1).astype(np.uint8)


def popcount(x):
    x = np.asarray(x, dtype=np.uint64)
    return (POPC[(x & np.uint64(0xFFFF)).astype(np.int64)].astype(np.int16)
            + POPC[((x >> np.uint64(16)) & np.uint64(0xFFFF)).astype(np.int64)]
            + POPC[((x >> np.uint64(32)) & np.uint64(0xFFFF)).astype(np.int64)]
            + POPC[((x >> np.uint64(48)) & np.uint64(0xFFFF)).astype(np.int64)])


# ============================================================ spectra
def spectrum(table, m, bits_per_res=None, want_top=0):
    """Exact Pauli-weight spectrum of a 2^m energy table.

    Returns a dict with Parseval / reconstruction checks, the weight shares V_w, the
    residue-order shares D_d (if bits_per_res given), and summary statistics.
    """
    t = np.asarray(table, float).ravel()
    assert t.size == 1 << m
    c = walsh_coeffs(t)
    idx = np.arange(1 << m, dtype=np.uint64)
    w = popcount(idx).astype(np.int64)
    c2 = c * c
    var = float(c2[1:].sum())
    out = {
        "m_bits": int(m),
        "mean": float(c[0]),
        "var": var,
        "std": float(np.sqrt(var)),
        "range": [float(t.min()), float(t.max())],
        # checks
        "parseval_abs_err": float(abs(float((c2).sum()) - float((t * t).mean()))),
        "parseval_rel_err": float(abs(float((c2).sum()) - float((t * t).mean()))
                                  / max(float((t * t).mean()), 1e-300)),
        "reconstruct_max_abs_err": float(np.abs(ifwht_reconstruct(c) - t).max()),
        "reconstruct_max_rel_err": float(np.abs(ifwht_reconstruct(c) - t).max()
                                         / max(float(np.abs(t).max()), 1e-300)),
    }
    if var <= 0.0:
        out["constant"] = True
        out["V"] = [0.0] * (m + 1)
        return out
    out["constant"] = False
    V = np.bincount(w[1:], weights=c2[1:], minlength=m + 1) / var
    out["V"] = [float(x) for x in V]
    out["W_mean"] = float((np.arange(m + 1) * V).sum())
    out["W_eff"] = float(1.0 / max(float((V ** 2).sum()), 1e-300))
    out["W_median"] = int(np.searchsorted(np.cumsum(V), 0.5))
    out["tail_gt1"] = float(V[2:].sum())
    out["tail_gt2"] = float(V[3:].sum())
    out["tail_gt3"] = float(V[4:].sum())
    out["V1"] = float(V[1])
    out["V2"] = float(V[2]) if m >= 2 else 0.0
    out["V3"] = float(V[3]) if m >= 3 else 0.0
    out["max_weight_with_mass"] = int(np.max(np.nonzero(V > 1e-14)[0])) if (V > 1e-14).any() else 0

    if bits_per_res:
        b = int(bits_per_res)
        L = m // b
        # residue-degree of each Walsh subset
        d = np.zeros(1 << m, dtype=np.int16)
        for r in range(L):
            mask = np.uint64(((1 << b) - 1) << (r * b))
            d += ((idx & mask) != 0).astype(np.int16)
        D = np.bincount(d[1:], weights=c2[1:], minlength=L + 1) / var
        out["bits_per_res"] = b
        out["n_res"] = int(L)
        out["D"] = [float(x) for x in D]
        out["D_mean"] = float((np.arange(L + 1) * D).sum())
        out["D_eff"] = float(1.0 / max(float((D ** 2).sum()), 1e-300))
        out["D1"] = float(D[1])
        out["D_tail_gt2"] = float(D[3:].sum())
        # exact order->weight bound check
        ok_lo = bool(np.all(w[1:] >= d[1:]))
        ok_hi = bool(np.all(w[1:] <= d[1:].astype(np.int64) * b))
        out["order_weight_bound_holds"] = bool(ok_lo and ok_hi)
        # mean weight conditional on residue order (the encoding's "spreading factor")
        num = np.bincount(d[1:], weights=c2[1:] * w[1:], minlength=L + 1)
        den = np.bincount(d[1:], weights=c2[1:], minlength=L + 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            out["W_given_D"] = [float(x) for x in np.where(den > 0, num / np.maximum(den, 1e-300),
                                                           np.nan)]
        out["spread"] = float(out["W_mean"] / out["D_mean"]) if out["D_mean"] > 0 else float("nan")
    if want_top:
        o = np.argsort(-c2[1:])[:want_top] + 1
        out["top_terms"] = [{"S": int(s), "weight": int(w[s]), "c": float(c[s]),
                             "share": float(c2[s] / var)} for s in o]
    return out


def check_order_weight_bound(m, b):
    """Independent re-derivation of the D->W bound, used as a unit test."""
    idx = np.arange(1 << m, dtype=np.uint64)
    w = popcount(idx).astype(np.int64)
    L = m // b
    d = np.zeros(1 << m, dtype=np.int64)
    for r in range(L):
        mask = np.uint64(((1 << b) - 1) << (r * b))
        d += ((idx & mask) != 0).astype(np.int64)
    return bool(np.all(w[1:] >= d[1:])) and bool(np.all(w[1:] <= d[1:] * b))


# ============================================================ encodings
def gray(i):
    return i ^ (i >> 1)


def perm_binary(k):
    return np.arange(k)


def perm_gray(k):
    """State j -> the bit pattern gray(j).  A relabelling of the k states."""
    return np.array([gray(j) for j in range(k)], dtype=int)


def apply_relabel(T, perms):
    """`T` is a (k,)*L table indexed by state indices; `perms[r][j]` is the BIT PATTERN
    assigned to state j of residue r.  Returns the table indexed by bit patterns."""
    out = T
    for r, p in enumerate(perms):
        inv = np.empty_like(np.asarray(p))
        inv[np.asarray(p)] = np.arange(len(p))
        out = np.take(out, inv, axis=r)
    return out


def onehot_table(T, k, fill="mean", penalty=None):
    """Embed a (k,)*L table into the 2^(L*k) one-hot register.

    `fill="mean"`  -> infeasible strings take the mean feasible energy (a neutral extension)
    `fill="penalty"` -> infeasible strings take mean + penalty * sum_r (popcount_r - 1)^2

    BOTH are modelling choices; the physical energy is defined ONLY on the feasible
    subspace, and the choice changes the spectrum.  That is itself the result.
    """
    L = T.ndim
    m = L * k
    assert m <= 22, "one-hot register too large to enumerate"
    base = float(np.asarray(T).mean())
    E = np.full(1 << m, base, dtype=float)
    idx = np.arange(1 << m, dtype=np.uint64)
    if fill == "penalty":
        viol = np.zeros(1 << m, dtype=np.float64)
        for r in range(L):
            mask = np.uint64(((1 << k) - 1) << (r * k))
            pc = popcount(idx & mask).astype(np.float64)
            viol += (pc - 1.0) ** 2
        E += float(penalty) * viol
    # write the feasible points: state j of residue r sets bit r*k + j
    flat = np.asarray(T).ravel()
    states = np.indices((k,) * L).reshape(L, -1).T              # (k^L, L), C order
    pos = np.zeros(len(states), dtype=np.int64)
    for r in range(L):
        pos += (np.int64(1) << (r * k + states[:, r]))
    E[pos] = flat
    return E


# ============================================================ sampling estimator
def krawtchouk_table(m):
    """K[w, a] = [t^w] (1+t)^{m-a} (1-t)^a  with a = number of DISAGREEING bits.

    Then  sum_{|S|=w} chi_S(x) chi_S(y) = K[w, hamming(x,y)].
    """
    K = np.zeros((m + 1, m + 1), dtype=np.float64)
    for a in range(m + 1):
        p = np.zeros(m + 1)
        p[0] = 1.0
        for _ in range(m - a):                # multiply by (1 + t)
            p[1:] += p[:-1].copy()
        for _ in range(a):                    # multiply by (1 - t)
            p[1:] -= p[:-1].copy()
        K[:, a] = p
    return K


def sampled_spectrum(X, f, m, blocks=8, chunk=256, seed=0):
    """UNBIASED estimator of sum_{|S|=w} c_S^2 from N sampled (bitstring, energy) pairs.

    X : (N,) uint64 bitstrings, drawn i.i.d. uniform on {0,1}^m.
    f : (N,) energies.

    Uses the U-statistic  V_w_hat = 1/(N(N-1)) sum_{a != b} f_a f_b K[w, ham(x_a,x_b)],
    which has expectation exactly sum_{|S|=w} c_S^2 for every w (including w=0, where it
    estimates the squared mean).  The pairwise sum is reduced to a Hamming histogram, so
    all m+1 weights come out of one O(N^2) pass.

    Returns the normalised V_w, the raw sums, and a block-split standard error.
    """
    X = np.asarray(X, dtype=np.uint64)
    f = np.asarray(f, dtype=float)
    N = len(X)
    K = krawtchouk_table(m)

    def hist_for(sub_idx):
        xs = X[sub_idx]
        fs = f[sub_idx]
        n = len(xs)
        H = np.zeros(m + 1, dtype=np.float64)
        for a0 in range(0, n, chunk):
            xa = xs[a0:a0 + chunk]
            fa = fs[a0:a0 + chunk]
            ham = popcount(xa[:, None] ^ xs[None, :]).astype(np.int64)
            wgt = fa[:, None] * fs[None, :]
            H += np.bincount(ham.ravel(), weights=wgt.ravel(), minlength=m + 1)[:m + 1]
        # remove the diagonal a == b
        H[0] -= float((fs * fs).sum())
        return H / (n * (n - 1))

    H = hist_for(np.arange(N))
    raw = K @ H                                   # raw[w] = sum_{|S|=w} c_S^2  (unbiased)
    rng = np.random.default_rng(seed)
    order = rng.permutation(N)
    parts = np.array_split(order, blocks)
    est = np.array([K @ hist_for(p) for p in parts])          # (blocks, m+1)
    var_tot = raw[1:].sum()
    out = {
        "N": int(N), "m_bits": int(m), "blocks": int(blocks),
        "raw": [float(x) for x in raw],
        "var_est": float(var_tot),
        "V": [float(x / var_tot) for x in raw],               # index 0 is the mean^2 share
    }
    Vb = est[:, 1:] / np.maximum(est[:, 1:].sum(1, keepdims=True), 1e-300)
    out["V_block_sd"] = [0.0] + [float(s) for s in Vb.std(0, ddof=1) / np.sqrt(blocks)]
    V = np.array(out["V"])
    V = np.concatenate([[0.0], V[1:] / max(V[1:].sum(), 1e-300)])
    out["V_norm"] = [float(x) for x in V]
    out["W_mean"] = float((np.arange(m + 1) * V).sum())
    out["V1"] = float(V[1]); out["V2"] = float(V[2]); out["V3"] = float(V[3])
    out["tail_gt3"] = float(V[4:].sum())
    Wb = (np.arange(m + 1)[None, :] * np.concatenate(
        [np.zeros((blocks, 1)), Vb], axis=1)).sum(1)
    out["W_mean_sd"] = float(Wb.std(ddof=1) / np.sqrt(blocks))
    return out


# ============================================================ energy tables
def enumerate_states(L, k):
    """All k^L configurations in the SAME order as the flattened (k,)*L table (C order)."""
    return np.indices((k,) * L).reshape(L, -1).T.astype(np.int64)


def space_for(pdb_id, L, k):
    """Truncate a tuning target to its first L residues (the geo/qarch convention)."""
    t = {x["pdb"]: x for x in Q.I.targets()}[pdb_id]
    return Q.Space(pdb_id, k, seq=t["seq"][:L], n=L, fold=t["fold"])


def legacy_tables(space, S, chunk=4096):
    """Per-term Legacy tables + the weighted total, for the given configurations."""
    comp = Q.legacy_components(space, S, chunk=chunk)
    tot = Q.legacy_total(comp)
    return comp, tot


def write(name, obj):
    return Q.write(name, obj)
