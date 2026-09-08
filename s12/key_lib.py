"""RETRIEVAL-KEY agent (sprint 12) -- shared library.

Question: can a SEQUENCE-PREDICTED local backbone conformation profile replace or augment
BLOSUM62 as the retrieval key?

Contents
    * torsion-bin alphabets (ABEGO-4, ABEGO-5, data-driven kmeans-K)
    * native torsions for the 126 tuning targets (from peptide_db)
    * key construction + retrieval (matched K=500, stable argsort like the shipped path)
    * per-fold LFO training corpora built from the universe files (leakage-safe by
      construction: universe(T) contains only out-of-fold peptides + this fold's fragments)
"""
from __future__ import annotations
import os, sys, json, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I

CACHE = I.CACHE
ALPHABET = I.ALPHABET

# --------------------------------------------------------------------- torsion alphabets
def abego4(phi, psi):
    """Standard Rosetta ABEGO minus the cis-omega 'O' state (omega unavailable).
    0=A (alpha)  1=B (beta/extended)  2=G (left-handed alpha)  3=E (left-handed extended)"""
    phi = np.asarray(phi, float); psi = np.asarray(psi, float)
    ph = np.degrees(np.arctan2(np.sin(phi), np.cos(phi)))
    ps = np.degrees(np.arctan2(np.sin(psi), np.cos(psi)))
    neg = ph < 0
    a = (ps > -75.0) & (ps <= 50.0)          # alpha psi band
    g = (ps > -100.0) & (ps <= 100.0)        # left-handed psi band
    out = np.where(neg, np.where(a, 0, 1), np.where(g, 2, 3))
    return out.astype(np.int8)


def abego5(phi, psi):
    """ABEGO with the B state split into P (polyproline-II, phi > -100) and B (beta)."""
    b = abego4(phi, psi).astype(np.int8)
    ph = np.degrees(np.arctan2(np.sin(np.asarray(phi, float)), np.cos(np.asarray(phi, float))))
    out = b.copy()
    out[(b == 1) & (ph > -100.0)] = 4        # 4 = P
    return out


_GRID_EDGES_PHI = np.array([-180.0, -120.0, -90.0, -60.0, 0.0, 180.0])
_GRID_EDGES_PSI = np.array([-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0])


def _wrap(x):
    x = np.asarray(x, float)
    return np.degrees(np.arctan2(np.sin(x), np.cos(x)))


class KMeansBins:
    """Data-driven Ramachandran binning: kmeans on (cos phi, sin phi, cos psi, sin psi).

    Fitted ONCE on a large, fold-independent sample of library torsions.  It is an
    unsupervised discretisation of the torsion space (no native / no pool involvement),
    so it is not a leakage channel; it is refit per fold anyway where cheap.
    """

    def __init__(self, k, seed=0):
        self.k = int(k); self.seed = int(seed); self.C = None

    @staticmethod
    def _feat(phi, psi):
        ph = np.asarray(phi, float).ravel(); ps = np.asarray(psi, float).ravel()
        return np.stack([np.cos(ph), np.sin(ph), np.cos(ps), np.sin(ps)], 1)

    def fit(self, phi, psi, iters=40):
        X = self._feat(phi, psi)
        rng = np.random.default_rng(self.seed)
        if len(X) > 200000:
            X = X[rng.permutation(len(X))[:200000]]
        C = X[rng.permutation(len(X))[: self.k]].copy()
        for _ in range(iters):
            d = ((X[:, None, :] - C[None]) ** 2).sum(2) if len(X) < 20000 else None
            if d is None:
                lab = np.empty(len(X), np.int32)
                for s in range(0, len(X), 20000):
                    blk = X[s:s + 20000]
                    lab[s:s + 20000] = np.argmin(((blk[:, None, :] - C[None]) ** 2).sum(2), 1)
            else:
                lab = np.argmin(d, 1)
            newC = C.copy()
            for c in range(self.k):
                m = lab == c
                if m.any():
                    newC[c] = X[m].mean(0)
            if np.allclose(newC, C, atol=1e-6):
                C = newC; break
            C = newC
        self.C = C
        return self

    def assign(self, phi, psi):
        X = self._feat(phi, psi)
        out = np.empty(len(X), np.int32)
        for s in range(0, len(X), 20000):
            blk = X[s:s + 20000]
            out[s:s + 20000] = np.argmin(((blk[:, None, :] - self.C[None]) ** 2).sum(2), 1)
        return out.reshape(np.asarray(phi).shape).astype(np.int8)


_KM_CACHE = {}


def kmeans_bins(k):
    """Global kmeans discretisation, fitted on a fixed sample of the library torsions."""
    if k in _KM_CACHE:
        return _KM_CACHE[k]
    path = os.path.join(CACHE, f"key_km{k}.npz")
    if os.path.exists(path):
        z = np.load(path); km = KMeansBins(k); km.C = z["C"]
        _KM_CACHE[k] = km; return km
    tg = I.targets()
    rng = np.random.default_rng(0)
    PH, PS = [], []
    for t in [tg[i] for i in rng.permutation(len(tg))[:12]]:
        u = I.load_univ(t["pdb"])
        sel = rng.permutation(len(u["PHI"]))[:4000]
        PH.append(u["PHI"][sel].ravel()); PS.append(u["PSI"][sel].ravel())
        del u
    km = KMeansBins(k).fit(np.concatenate(PH), np.concatenate(PS))
    np.savez(path, C=km.C)
    _KM_CACHE[k] = km
    return km


def ss3(phi, psi):
    """3-state H/E/C from torsions (the coordinator's resolution control -- this is the
    alphabet the forensics agent found nearly worthless as a key).
    0=H 1=E 2=C."""
    ph = _wrap(phi); ps = _wrap(psi)
    H = (ph > -160.0) & (ph < -20.0) & (ps > -120.0) & (ps < 50.0) & \
        (ph > -100.0) & (ps > -80.0) & (ps < 5.0)
    H = (ph > -100.0) & (ph < -30.0) & (ps > -80.0) & (ps < -5.0)
    E = (ph > -180.0) & (ph < -45.0) & (((ps > 90.0) & (ps <= 180.0)) | (ps < -170.0))
    out = np.full(np.shape(ph), 2, np.int8)
    out[E] = 1
    out[H] = 0
    return out


ALPHABETS = {
    "ss3": (3, ss3),
    "abego4": (4, abego4),
    "abego5": (5, abego5),
    "km8": (8, lambda p, s: kmeans_bins(8).assign(p, s)),
    "km16": (16, lambda p, s: kmeans_bins(16).assign(p, s)),
    "km12": (12, lambda p, s: kmeans_bins(12).assign(p, s)),
}


def nbins(name):
    return ALPHABETS[name][0]


def bins_of(name, phi, psi):
    return ALPHABETS[name][1](phi, psi)


# --------------------------------------------------------------------- native torsions
_NAT = None


def native_torsions():
    """{pdb: (phi (n,), psi (n,))} for every tuning target, from peptide_db (ORACLE)."""
    global _NAT
    if _NAT is not None:
        return _NAT
    path = os.path.join(CACHE, "key_native_torsions.npz")
    if os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        _NAT = {str(k): (np.asarray(v[0], float), np.asarray(v[1], float))
                for k, v in zip(z["pdbs"], z["tor"])}
        return _NAT
    import peptide_db as P
    db = {p.pdb.upper(): p for p in P.load()}
    out = {}
    for t in I.targets():
        p = db.get(t["pdb"].upper())
        if p is None or p.seq != t["seq"]:
            # fall back to sequence match
            cand = [q for q in db.values() if q.seq == t["seq"]]
            p = cand[0] if cand else None
        if p is None:
            continue
        out[t["pdb"]] = (np.asarray(p.phi, float), np.asarray(p.psi, float))
    np.savez(path, pdbs=np.array(list(out), object),
             tor=np.array([[v[0], v[1]] for v in out.values()], object))
    _NAT = out
    return out


# --------------------------------------------------------------------- retrieval
def zscore(x):
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / (s if s > 1e-9 else 1.0)


def retrieve(score, k=I.K):
    """Stable top-k by DESCENDING score, matching the shipped `argsort(-sim)` convention."""
    order = np.argsort(-np.asarray(score, float), kind="stable")
    return np.asarray(order[:k], int)


def hard_key_score(win_bins, tgt_bins):
    """Number of matching bins (the literature agent's key)."""
    return (win_bins == np.asarray(tgt_bins)[None, :]).sum(1).astype(float)


def soft_key_score(win_bins, post):
    """sum_i log P(bin_i = win_bins[:,i]); post is (n, nb) row-stochastic."""
    lp = np.log(np.clip(np.asarray(post, float), 1e-6, 1.0))
    n = win_bins.shape[1]
    return lp[np.arange(n)[None, :], win_bins].sum(1)


# --------------------------------------------------------------------- evaluation
def pool_stats(rr, idx, band_ref=None):
    r = np.asarray(rr, float)[np.asarray(idx, int)]
    ref = float(np.asarray(rr, float).min()) if band_ref is None else float(band_ref)
    return {"best": float(r.min()), "mean": float(r.mean()),
            "band": int((r <= ref + I.BAND).sum()), "p10": float(np.percentile(r, 10))}


def emit(u, idx, dg, m=I.M, lam=0.3):
    """Full production chain on a pool: shipped-score top-m -> coordinate average -> project.
    Returns dict with rmsd at lam=0 (fit) and lam (projected), plus top-m diagnostics."""
    idx = np.asarray(idx, int)
    i, j = I.pair_index(u["n"])
    D = I.pair_dists(u["W"][idx], i, j)
    sc = I.shipped_score(dg, D.astype(np.float32).astype(float))
    sub = np.argsort(sc, kind="stable")[:m]
    W = u["W"][idx[sub]]
    C, _ = I.coordinate_average(W)
    out = I.project(C, u["seq"], u["fold"], lam=lam)
    rr = u["rr"][idx]
    return {"fit": I.ca_rmsd(out["fit_ca"], u["nat_ca"]),
            "proj": I.ca_rmsd(out["ca"], u["nat_ca"]),
            "argmin": float(rr[int(np.argmin(sc))]),
            "top_m_best": float(rr[sub].min()),
            "top_m_mean": float(rr[sub].mean()),
            "band_recall_in_m": int((rr[sub] <= rr.min() + I.BAND).sum())}
