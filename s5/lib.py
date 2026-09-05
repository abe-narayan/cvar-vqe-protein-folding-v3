"""Shared primitives for the Sprint 5 retrieval experiments.

Factored out of `s5/retrieval_probe.py` because importing that module to reuse its
helpers re-ran the entire probe as a side effect -- it had no __main__ guard, so the
first consumer silently spent minutes recomputing a finished experiment.
Everything here is pure: no experiment runs on import.
"""
import json
import os
import sys

sys.path.insert(0, os.path.abspath("."))
import numpy as np

import distogram as dgm
import peptide_db as db

DEV = ("1AL1,1GJE,1KVG,1MPV,1YY1,2DCX,2DX2,2FXZ,2MI1,2MOC,2N08,2NS4,"
       "5XER,5ZGD,6GS5,6M19,6RSM,7JGY,7K1M,7M77,7YSS,8DHZ,9U4Q,9U7M").split(",")
KS = (1, 5, 20, 100, 500, 2000)
OUT = "s5/retrieval_probe.json"

AA = "ARNDCQEGHILKMFPSTWYV"
#: BLOSUM62, upper triangle by AA order above, flattened.
_B62 = """4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
4 2 -3 1 0 -3 -2 -1 -3 -1 3
4 -2 2 0 -3 -2 -1 -2 -1 1
5 -1 -3 -1 0 -1 -3 -2 -2
5 0 -2 -1 -1 -1 -1 1
6 -4 -2 -2 1 3 -1
7 -1 -1 -4 -3 -2
4 1 -3 -2 -2
5 -2 -2 0
11 2 -3
7 -1
4"""


def _blosum():
    M = np.zeros((20, 20))
    rows = [r.split() for r in _B62.strip().split("\n")]
    for i, r in enumerate(rows):
        for k, v in enumerate(r):
            j = i + k
            M[i, j] = M[j, i] = float(v)
    return M


B62 = _blosum()
IDX = {a: i for i, a in enumerate(AA)}


def encode(seq):
    return np.array([IDX.get(c, 0) for c in seq], dtype=np.int64)


def windows(pool, n):
    """(W,n,3) CA windows and (W,n) encoded sequences."""
    cas, seqs = [], []
    for q in pool:
        m = len(q.seq)
        if m < n:
            continue
        e = encode(q.seq)
        for s in range(m - n + 1):
            cas.append(q.ca[s:s + n])
            seqs.append(e[s:s + n])
    return np.stack(cas), np.stack(seqs)


def kabsch_rmsd_batch(W, T):
    W = W - W.mean(1, keepdims=True)
    T = T - T.mean(0, keepdims=True)
    H = np.einsum("bni,nj->bij", W, T)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("bij,bjk->bik", Vt.transpose(0, 2, 1),
                                        U.transpose(0, 2, 1))))
    S = S.copy(); S[:, -1] *= d
    num = (W ** 2).sum((1, 2)) + (T ** 2).sum() - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(num, 0.0) / W.shape[1])



KS = (1, 5, 20, 100, 500, 2000)
OUT = "s5/retrieval_probe.json"


def windows_full(pool, n):
    """(W,n,3) CA windows, (W,n) encoded sequence, (W,n) phi, (W,n) psi, and provenance.

    A window carries its parent's torsions directly, so a retrieved candidate has a real
    conformation with no reconstruction step and nothing to infer. phi[s] of a window
    starting at s is defined by the residue before it in the parent chain, which exists;
    the window is a genuine fragment conformation, not a truncation artefact.
    """
    cas, seqs, phis, psis, src = [], [], [], [], []
    for k, q in enumerate(pool):
        m = len(q.seq)
        if m < n:
            continue
        e = encode(q.seq)
        for s in range(m - n + 1):
            cas.append(q.ca[s:s + n]); seqs.append(e[s:s + n])
            phis.append(q.phi[s:s + n]); psis.append(q.psi[s:s + n])
            src.append((k, s))
    return (np.stack(cas), np.stack(seqs), np.stack(phis).astype(np.float32),
            np.stack(psis).astype(np.float32), src)
