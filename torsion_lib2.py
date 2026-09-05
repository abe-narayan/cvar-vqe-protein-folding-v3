"""Per-residue torsion state libraries built from the peptide database.

The shipped encoding gives every residue one of four *class* libraries (GENERAL, GLY,
PRO, PRE_PRO). That spends its k states on the union of everything a class ever does. A
residue's actual (phi, psi) distribution is much tighter than its class's once you
condition on the residue type and its neighbours, and the states are free to differ per
residue -- the bitstring layout only cares that every residue offers exactly k of them.

So the library here is *per residue*: for residue i of a target sequence, the states are
k-means centres over the (phi, psi) pairs observed in the database for residues matching
i's context, with a back-off chain from the most specific context to the class pool. Same
qubit count, lower representation floor.

Contexts, most specific first:
    (aa[i-1], aa[i], aa[i+1])   exact triplet
    (aa[i], aa[i+1])            and (aa[i-1], aa[i])
    (aa[i],)                    residue type
    class                       GENERAL / GLY / PRO / PRE_PRO

A context is used when it has at least `MIN_OBS` observations, and observations from
weaker contexts are always blended in with a smaller weight so a specific context with
just enough data cannot produce k states that are all one cluster.
"""
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import peptide_db as pdb
import representations as reps

_CLASSES = (reps.CLASS_GENERAL, reps.CLASS_GLY, reps.CLASS_PRO,
            reps.CLASS_PRE_PRO)
MIN_OBS = 40
#: Weight of the backed-off pool relative to the specific one. Not zero: with 60
#: observations, k=16 clusters would otherwise be fitted to ~4 points each.
BACKOFF_W = 0.25


def _circ_kmeans(X4: np.ndarray, w: np.ndarray, k: int, seed: int = 1,
                 iters: int = 60) -> np.ndarray:
    """Weighted k-means on the 4-D (cos phi, sin phi, cos psi, sin psi) embedding."""
    m = len(X4)
    k = min(k, m)
    rng = np.random.default_rng(seed)
    # k-means++ on the weights, so rare-but-real basins are not missed by a uniform draw.
    idx = [int(rng.choice(m, p=w / w.sum()))]
    d2 = ((X4 - X4[idx[0]]) ** 2).sum(1)
    for _ in range(k - 1):
        p = d2 * w
        s = p.sum()
        idx.append(int(rng.choice(m, p=p / s)) if s > 0 else int(rng.integers(m)))
        d2 = np.minimum(d2, ((X4 - X4[idx[-1]]) ** 2).sum(1))
    C = X4[idx].copy()
    for _ in range(iters):
        lab = ((X4[:, None, :] - C[None]) ** 2).sum(-1).argmin(1)
        newC = C.copy()
        for j in range(k):
            sel = lab == j
            if sel.any():
                ww = w[sel][:, None]
                newC[j] = (X4[sel] * ww).sum(0) / ww.sum()
        if np.allclose(newC, C):
            break
        C = newC
    return C


def _angles_from_centres(C: np.ndarray) -> np.ndarray:
    return np.column_stack([np.arctan2(C[:, 1], C[:, 0]),
                            np.arctan2(C[:, 3], C[:, 2])])


class ContextLibrary:
    """Observation pools keyed by sequence context, built once per holdout database."""

    def __init__(self, entries: Sequence[pdb.Peptide]):
        self.pools: Dict[tuple, List[Tuple[float, float]]] = {}
        for p in entries:
            cls = reps.residue_classes(p.seq, p.n)
            for i in range(p.n):
                a = p.seq[i]
                l = p.seq[i - 1] if i > 0 else "^"
                r = p.seq[i + 1] if i + 1 < p.n else "$"
                t = (float(p.phi[i]), float(p.psi[i]))
                for key in ((l, a, r), (a, r), (l, a), (a,), (cls[i],)):
                    self.pools.setdefault(key, []).append(t)
        self._cache: Dict[tuple, np.ndarray] = {}

    def _keys(self, seq: str, i: int, cls: str) -> List[tuple]:
        a = seq[i]
        l = seq[i - 1] if i > 0 else "^"
        r = seq[i + 1] if i + 1 < len(seq) else "$"
        return [(l, a, r), (a, r), (l, a), (a,), (cls,)]

    def states(self, seq: str, i: int, cls: str, k: int, seed: int = 1) -> np.ndarray:
        """(k, 2) array of (phi, psi) in radians for residue `i` of `seq`."""
        key = (seq[max(0, i - 1):i + 2], seq[i], cls, k, i == 0, i == len(seq) - 1)
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        chain = self._keys(seq, i, cls)
        pts: List[Tuple[float, float]] = []
        wts: List[float] = []
        used = False
        for depth, kk in enumerate(chain):
            pool = self.pools.get(kk)
            if not pool:
                continue
            if not used and len(pool) >= MIN_OBS:
                pts.extend(pool)
                wts.extend([1.0] * len(pool))
                used = True
            elif used:
                pts.extend(pool)
                wts.extend([BACKOFF_W] * len(pool))
                break
        if not used:                      # fall back to the class pool outright
            pool = self.pools.get((cls,)) or self.pools.get(("GENERAL",)) or []
            pts, wts = list(pool), [1.0] * len(pool)
        P = np.asarray(pts, float)
        w = np.asarray(wts, float)
        X4 = np.column_stack([np.cos(P[:, 0]), np.sin(P[:, 0]),
                              np.cos(P[:, 1]), np.sin(P[:, 1])])
        ang = _angles_from_centres(_circ_kmeans(X4, w, k, seed=seed))
        if len(ang) < k:                  # pad by repeating (degenerate tiny pool)
            ang = np.vstack([ang, np.repeat(ang[:1], k - len(ang), axis=0)])
        self._cache[key] = ang
        return ang


@lru_cache(maxsize=32)
def _class_library(k: int, exclude_seq: str, seed: int) -> tuple:
    """k states per residue CLASS, clustered over the held-out database."""
    entries = pdb.holdout(exclude_seq) if exclude_seq else pdb.load()
    lib = ContextLibrary(entries)
    out = {}
    for c in _CLASSES:
        pool = np.asarray(lib.pools.get((c,))
                          or lib.pools[(reps.CLASS_GENERAL,)], float)
        X4 = np.column_stack([np.cos(pool[:, 0]), np.sin(pool[:, 0]),
                              np.cos(pool[:, 1]), np.sin(pool[:, 1])])
        out[c] = _angles_from_centres(
            _circ_kmeans(X4, np.ones(len(X4)), k, seed=seed))
    return tuple(sorted(out.items()))


@lru_cache(maxsize=32)
def library_for(sequence: str, k: int, exclude_seq: str = "", seed: int = 1,
                mode: str = "class") -> np.ndarray:
    """``(n, k, 2)`` per-residue (phi, psi) table in radians, holding out `exclude_seq`.

    ``mode="class"`` gives every residue the k states of its class (GENERAL/GLY/PRO/
    PRE_PRO); ``mode="ctx"`` gives each residue its own states from the back-off context
    chain above. The context form was measured against the class form on 35 targets at
    k = 4, 8, 16, 32 and is NOT better -- mean floor 1.39/1.02/0.76/0.58 A against
    1.38/0.93/0.69/0.54. Conditioning on the triplet fragments the observation pool faster
    than it sharpens the distribution, so the class form is the default and the context
    form is kept because the comparison is the evidence for that choice.
    """
    cls = reps.residue_classes(sequence, len(sequence))
    if mode == "class":
        lib = dict(_class_library(k, exclude_seq, seed))
        return np.stack([lib[c] for c in cls])
    entries = pdb.holdout(exclude_seq) if exclude_seq else pdb.load()
    ctx = ContextLibrary(entries)
    return np.stack([ctx.states(sequence, i, cls[i], k, seed=seed)
                     for i in range(len(sequence))])


class PerResidueTorsion(reps.TorsionStateRepresentation):
    """`TorsionStateRepresentation` with an explicit per-residue state table.

    Everything downstream -- bitstring layout, chi1 bits, `build_all`, the Hamiltonians,
    `batch_features` -- reads `_phi` / `_psi`, so supplying the table is all it takes.
    """

    name = "torsion_ctx"

    def __init__(self, sequence: str, table: np.ndarray, chi_bits: bool = True):
        n, k, _ = table.shape
        if n != len(sequence):
            raise ValueError(f"table has {n} residues, sequence has {len(sequence)}")
        reps.STATE_LIBRARIES.setdefault(
            k, {c: [(-63.0, -42.0)] * k for c in
                (reps.CLASS_GENERAL, reps.CLASS_GLY, reps.CLASS_PRO,
                 reps.CLASS_PRE_PRO)})
        super().__init__(n, n_states=k, sequence=sequence, chi_bits=chi_bits)
        self._phi = np.ascontiguousarray(table[:, :, 0])
        self._psi = np.ascontiguousarray(table[:, :, 1])
        self.libraries = [[(float(np.degrees(a)), float(np.degrees(b)))
                           for a, b in table[i]] for i in range(n)]


def make(sequence: str, k: int = 8, exclude_seq: Optional[str] = None,
         chi_bits: bool = False, seed: int = 1,
         mode: str = "class") -> PerResidueTorsion:
    tab = library_for(sequence, k, exclude_seq or "", seed, mode)
    return PerResidueTorsion(sequence, tab, chi_bits=chi_bits)
