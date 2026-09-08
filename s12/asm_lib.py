"""Sprint 12 -- ASSEMBLY agent shared machinery.

Sub-length window banks over EXACTLY the legal library of a fold (out-of-fold peptides +
`core.predict._fold_fragments(fold, 5)`), a fast many-vs-many Kabsch, and a batched
ideal-geometry CA builder wrapper.

Everything here is leakage-safe by construction: a bank is keyed by (fold, L) and the
peptide list is filtered with `peptide_db.folds(5)[seq] != fold`, which already removes the
target itself (its own sequence is in its own fold).
"""
from __future__ import annotations
import os, sys, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
BANKS = os.path.join(ROOT, "s12", "cache", "asm_banks")
os.makedirs(BANKS, exist_ok=True)

ALPHABET = "ARNDCQEGHILKMFPSTWYV"
_CODE = {c: k for k, c in enumerate(ALPHABET)}


def encode(seq):
    from s7 import audit
    return audit.encode(seq)


def B62():
    from s7 import audit
    return audit.B62


# ------------------------------------------------------------------ banks
_MEM = {}


def _legal_library(fold):
    """(entries, is_peptide) -- exactly the library `s8/generate.py:stage_univ` uses."""
    import peptide_db as db
    from core import predict as dgm
    folds = db.folds(5)
    peps = [q for q in db.load() if folds[q.seq] != fold]
    fr = list(dgm._fold_fragments(fold, 5))
    return peps + fr, [True] * len(peps) + [False] * len(fr)


def bank(fold, L, verbose=False):
    """All length-L windows of fold `fold`'s legal library.

    dict: W (nw,L,3) f64 CA, PHI/PSI (nw,L) f64, S (nw,L) int8, org (nw,) bool,
          src (nw,) int32 parent index, pos (nw,) int16 window start in parent.
    """
    key = (int(fold), int(L))
    if key in _MEM:
        return _MEM[key]
    path = os.path.join(BANKS, f"f{fold}_L{L}.npz")
    if os.path.exists(path):
        z = np.load(path)
        b = {k: z[k] for k in z.files}
    else:
        ents, isp = _legal_library(fold)
        cas, phis, psis, seqs, org, src, pos = [], [], [], [], [], [], []
        for k, (q, flag) in enumerate(zip(ents, isp)):
            m = len(q.seq)
            if m < L:
                continue
            e = encode(q.seq)
            for s in range(m - L + 1):
                cas.append(q.ca[s:s + L]); phis.append(q.phi[s:s + L])
                psis.append(q.psi[s:s + L]); seqs.append(e[s:s + L])
                org.append(flag); src.append(k); pos.append(s)
        b = dict(W=np.stack(cas).astype(np.float32), PHI=np.stack(phis).astype(np.float32),
                 PSI=np.stack(psis).astype(np.float32), S=np.stack(seqs).astype(np.int8),
                 org=np.array(org, bool), src=np.array(src, np.int32), pos=np.array(pos, np.int16))
        np.savez_compressed(path, **b)
        if verbose:
            print(f"  bank f{fold} L{L}: {len(b['W'])} windows", flush=True)
    b = {k: (np.asarray(v, np.float64) if k in ("W", "PHI", "PSI") else v) for k, v in b.items()}
    b["Wc"] = b["W"] - b["W"].mean(1, keepdims=True)
    b["n2"] = (b["Wc"] ** 2).sum((1, 2))
    _MEM[key] = b
    return b


def bank_ideal(fold, L, chunk=20000, lean=False):
    """Ideal-geometry CA rebuilt from each bank window's own real (phi,psi).

    Adds Ic (centred) and In2 to the bank dict; cached on disk.  This is what a window
    contributes to a CHAIN-CONSISTENT assembly, and its RMSD to a native segment is an
    exact lower bound on any assembly containing that piece in that slot.
    """
    b = bank(fold, L)
    if "Ic" in b:
        return b
    path = os.path.join(BANKS, f"ideal_f{fold}_L{L}.npz")
    if os.path.exists(path):
        Ideal = np.asarray(np.load(path)["I"], np.float64)
    else:
        parts = []
        for a in range(0, len(b["PHI"]), chunk):
            parts.append(build_ca(b["PHI"][a:a + chunk], b["PSI"][a:a + chunk]))
        Ideal = np.concatenate(parts, 0)
        np.savez_compressed(path, I=Ideal.astype(np.float32))
    b["Ic"] = Ideal - Ideal.mean(1, keepdims=True)
    b["In2"] = (b["Ic"] ** 2).sum((1, 2))
    if lean:
        b["nw"] = len(b["Ic"])
        for k in ("W", "Wc", "n2"):
            b.pop(k, None)
    return b


def drop_bank(fold=None, L=None):
    for k in list(_MEM):
        if (fold is None or k[0] == fold) and (L is None or k[1] == L):
            del _MEM[k]


# ------------------------------------------------------------------ fast Kabsch
def _sym3_eigvals(M):
    """Ascending-free analytic eigenvalues (e1>=e2>=e3) of a batch of symmetric 3x3."""
    m00 = M[..., 0, 0]; m11 = M[..., 1, 1]; m22 = M[..., 2, 2]
    m01 = M[..., 0, 1]; m02 = M[..., 0, 2]; m12 = M[..., 1, 2]
    q = (m00 + m11 + m22) / 3.0
    p1 = m01 ** 2 + m02 ** 2 + m12 ** 2
    p2 = (m00 - q) ** 2 + (m11 - q) ** 2 + (m22 - q) ** 2 + 2.0 * p1
    p = np.sqrt(np.maximum(p2, 0.0) / 6.0)
    pp = np.where(p > 1e-300, p, 1.0)
    b00 = (m00 - q) / pp; b11 = (m11 - q) / pp; b22 = (m22 - q) / pp
    b01 = m01 / pp; b02 = m02 / pp; b12 = m12 / pp
    det = (b00 * (b11 * b22 - b12 * b12) - b01 * (b01 * b22 - b12 * b02)
           + b02 * (b01 * b12 - b11 * b02))
    r = np.clip(det / 2.0, -1.0, 1.0)
    ph = np.arccos(r) / 3.0
    e1 = q + 2.0 * p * np.cos(ph)
    e3 = q + 2.0 * p * np.cos(ph + 2.0 * math.pi / 3.0)
    e2 = 3.0 * q - e1 - e3
    flat = p <= 1e-300
    if flat.any():
        e1 = np.where(flat, q, e1); e2 = np.where(flat, q, e2); e3 = np.where(flat, q, e3)
    return e1, e2, e3


def msd_many(Wc, n2, T, chunk=200000):
    """Mean-square CA deviation of every centred window in Wc (nw,L,3) against every
    query in T (m,L,3), after optimal superposition.  Returns (nw, m) float64.

    Analytic 3x3 singular values -- ~8x faster than np.linalg.svd and agrees to <1e-9 A.
    """
    Wc = np.asarray(Wc, float); T = np.asarray(T, float)
    if T.ndim == 2:
        T = T[None]
    L = Wc.shape[1]
    Tc = T - T.mean(1, keepdims=True)
    t2 = (Tc ** 2).sum((1, 2))
    out = np.empty((len(Wc), len(Tc)))
    for a in range(0, len(Wc), chunk):
        w = Wc[a:a + chunk]
        H = np.einsum("bli,qlj->bqij", w, Tc, optimize=True)
        M = np.einsum("bqki,bqkj->bqij", H, H, optimize=True)
        e1, e2, e3 = _sym3_eigvals(M)
        det = (H[..., 0, 0] * (H[..., 1, 1] * H[..., 2, 2] - H[..., 1, 2] * H[..., 2, 1])
               - H[..., 0, 1] * (H[..., 1, 0] * H[..., 2, 2] - H[..., 1, 2] * H[..., 2, 0])
               + H[..., 0, 2] * (H[..., 1, 0] * H[..., 2, 1] - H[..., 1, 1] * H[..., 2, 0]))
        s = (np.sqrt(np.maximum(e1, 0.0)) + np.sqrt(np.maximum(e2, 0.0))
             + np.sign(det) * np.sqrt(np.maximum(e3, 0.0)))
        out[a:a + chunk] = np.maximum(n2[a:a + chunk, None] + t2[None, :] - 2.0 * s, 0.0) / L
    return out


def rmsd_many(Wc, n2, T, chunk=200000):
    return np.sqrt(msd_many(Wc, n2, T, chunk))


# ------------------------------------------------------------------ CA builder
def build_ca(phi, psi):
    """Batched ideal-geometry CA trace; identical to `s12.instrument.build_ca`."""
    from core import project as pj
    phi = np.atleast_2d(np.asarray(phi, float)); psi = np.atleast_2d(np.asarray(psi, float))
    return np.asarray(pj.build_ca_exact(phi, psi), float)


def rmsd_to(X, T):
    """(b,n,3) vs (n,3) -> (b,) CA-RMSD."""
    X = np.asarray(X, float)
    Xc = X - X.mean(1, keepdims=True)
    n2 = (Xc ** 2).sum((1, 2))
    return np.sqrt(msd_many(Xc, n2, np.asarray(T, float)[None])[:, 0])
