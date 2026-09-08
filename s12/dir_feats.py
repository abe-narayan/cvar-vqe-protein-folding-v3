"""D3b. Deployable per-pair features for the error-SIGN head.

Everything here is computable at inference time from the target sequence and the shipped
leave-fold-out distogram.  No native coordinates, no `rr`, no pool statistics.
"""
from __future__ import annotations
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

_PEP = None
_HOT = None


def esm_for(seq):
    """(pca32 (n,32), con (n,n)) from the compact peptide bank / hot cache."""
    global _PEP, _HOT
    if _PEP is None:
        z = np.load(os.path.join(ROOT, "s12", "cache", "dir_esm_pep.npz"), allow_pickle=True)
        _PEP = {str(s): (np.asarray(a, np.float32), np.asarray(c, np.float32))
                for s, a, c in zip(z["seqs"], z["pca32"], z["con"])}
    if seq in _PEP:
        return _PEP[seq]
    from core import data as cdata
    return np.asarray(cdata.esm_embed(seq), np.float32), np.asarray(cdata.esm_contacts(seq), np.float32)


def _shell_resid(v, sep):
    out = np.zeros_like(v)
    for s in np.unique(sep):
        m = sep == s
        out[m] = v[m] - v[m].mean()
    return out


def _row_mean(v, i, j, n):
    """Mean of a per-pair quantity over every pair that touches each residue."""
    acc = np.zeros(n); cnt = np.zeros(n)
    np.add.at(acc, i, v); np.add.at(cnt, i, 1.0)
    np.add.at(acc, j, v); np.add.at(cnt, j, 1.0)
    return acc / np.maximum(cnt, 1.0)


NAMES = None


def build(seq, prob, exp, sd, centres):
    """X (npairs, D) of deployable features.  `prob` (npairs,17), `exp`/`sd` (npairs,)."""
    global NAMES
    from core import predict as cp
    n = len(seq)
    i, j = np.triu_indices(n, k=2)
    sep = (j - i).astype(float)
    prob = np.asarray(prob, np.float64); prob = prob / np.maximum(prob.sum(1, keepdims=True), 1e-12)
    exp = np.asarray(exp, np.float64); sd = np.asarray(sd, np.float64)
    C = np.asarray(centres, np.float64)

    # ---- shape of the predicted distribution
    ent = -(prob * np.log(prob + 1e-9)).sum(1)
    amax = prob.argmax(1)
    pmax = prob.max(1)
    mode = C[amax]
    cdf = np.cumsum(prob, 1)
    def q(p):
        k = np.clip((cdf < p).sum(1), 0, len(C) - 1)
        return C[k]
    q25, q50, q75 = q(0.25), q(0.5), q(0.75)
    skew = (exp - q50) / np.maximum(sd, 1e-3)
    # second mode: mass outside a +-2 bin window round the argmax
    win = np.abs(np.arange(len(C))[None, :] - amax[:, None]) <= 2
    tail = (prob * (~win)).sum(1)
    p_lt = np.column_stack([(prob * (C[None, :] < t)).sum(1) for t in (6.0, 8.0, 10.0, 12.0)])

    # ---- within-target context
    sr_exp = _shell_resid(exp, sep); sr_sd = _shell_resid(sd, sep)
    rexp_i = _row_mean(exp, i, j, n)[i]; rexp_j = _row_mean(exp, i, j, n)[j]
    rsd_i = _row_mean(sd, i, j, n)[i]; rsd_j = _row_mean(sd, i, j, n)[j]
    # implied radius of gyration from the predicted matrix (short pairs at ideal geometry)
    Dm = np.zeros((n, n)); Dm[i, j] = exp; Dm[j, i] = exp
    k1 = np.arange(n - 1); Dm[k1, k1 + 1] = Dm[k1 + 1, k1] = 3.80
    rg = float(np.sqrt(max((Dm ** 2).sum() / (2.0 * n * n), 0.0)))
    rg_ref = 2.2 * n ** 0.38                     # crude expectation for a peptide of length n
    gexp = float(exp.mean()); gsd = float(sd.mean()); gent = float(ent.mean())

    # ---- ESM
    E, con = esm_for(seq)
    con = np.asarray(con, np.float64)
    deg = con.sum(1) / max(1, n - 1)
    s_ax = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    shells = np.stack([np.where((s_ax >= a) & (s_ax < b), con, 0.0).sum(1) / max(1, n)
                       for a, b in ((3, 5), (5, 9), (9, 100))], axis=1)
    nb = np.stack([con[np.clip(i + a, 0, n - 1), np.clip(j + b, 0, n - 1)]
                   for a, b in ((-1, -1), (1, 1), (-1, 1), (1, -1))], axis=1)

    P = np.asarray(cp.residue_props(seq), np.float64)

    cols, names = [], []
    def add(v, nm):
        v = np.asarray(v, np.float64)
        if v.ndim == 1:
            v = v[:, None]; nms = [nm]
        else:
            nms = [f"{nm}{k}" for k in range(v.shape[1])]
        cols.append(v); names.extend(nms)

    add(sep, "sep"); add(np.log(sep), "logsep"); add(sep / n, "sepfrac")
    add(np.full(sep.shape, float(n)), "n")
    add(i / n, "ifrac"); add(j / n, "jfrac")
    add(np.minimum(i, n - 1 - j) / n, "endfrac")
    add(exp, "exp"); add(sd, "sd"); add(sd / np.maximum(exp, 1e-3), "sdrel")
    add(exp / np.maximum(sep, 1.0), "exp_per_sep")
    add(ent, "ent"); add(pmax, "pmax"); add(mode, "mode"); add(exp - mode, "exp_m_mode")
    add(q25, "q25"); add(q50, "q50"); add(q75, "q75"); add(q75 - q25, "iqr")
    add(skew, "skew"); add(tail, "tailmass"); add(p_lt, "plt")
    add(sr_exp, "shres_exp"); add(sr_sd, "shres_sd")
    add(rexp_i, "rowexp_i"); add(rexp_j, "rowexp_j")
    add(rsd_i, "rowsd_i"); add(rsd_j, "rowsd_j")
    add(np.full(sep.shape, rg), "rg"); add(np.full(sep.shape, rg / rg_ref), "rg_rel")
    add(np.full(sep.shape, gexp), "gexp"); add(np.full(sep.shape, gsd), "gsd")
    add(np.full(sep.shape, gent), "gent")
    add(exp - gexp, "exp_m_gexp")
    add(con[i, j], "con"); add(nb, "connb")
    add(deg[i], "deg_i"); add(deg[j], "deg_j")
    add(shells[i], "sh_i"); add(shells[j], "sh_j")
    add(P[i], "pi"); add(P[j], "pj")
    add(E[i], "ei"); add(E[j], "ej")
    X = np.hstack(cols).astype(np.float32)
    NAMES = names
    return X, names


SEP_COLS = ("sep", "logsep", "sepfrac", "n")
