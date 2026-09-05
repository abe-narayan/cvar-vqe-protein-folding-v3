"""CA-trace geometry prior: the term that tells a structure from its mirror image.

**The defect this fixes.** The search objective is dominated by a predicted CA-CA distance
matrix, and *a distance matrix is invariant under reflection*. The mirror image of any
candidate has an identical distance matrix and therefore an identical score. This is not a
tuning problem, it is a blind spot in the functional form, and it is expensive: a
per-target diagnosis found that 15-24% of the candidate pool for 1LE0 consists of
left-handed copies, and on two of three seeds the pipeline **selected one** (picked 3.32 A
where its mirror was 2.68 A; 3.54 A where the mirror was 2.24 A). 1K43 shows the same at
seed 2 (2.79 A picked, 1.34 A mirror).

**The fix.** Score the CA trace on the two internal coordinates that describe it: the
pseudo-bond-angle `theta_i` over (CA_i, CA_i+1, CA_i+2), which is mirror-invariant and
constrains local compactness, and the pseudo-torsion `tau_i` over four consecutive CA
atoms, which **changes sign under reflection**. Their joint distribution in deposited
structures is sharply bimodal -- right-handed helix near (89 deg, +50 deg), extended near
(120 deg, +/-180 deg) -- and the left-handed mirror of a helix lands where almost no real
protein sits. A joint histogram over (theta, tau), estimated from the held-out database and
consumed as `-log P`, therefore rejects mirrors on physics rather than on a fitted
constant.

This is a coarse-grained, residue-level, non-all-atom potential: it belongs to the Legacy
energy track and is kept conceptually distinct from the all-atom Amber track.

Estimated once per holdout and cached; scoring a batch is a digitise plus a gather.
"""
import os
from functools import lru_cache
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

import peptide_db as pdb

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "catrace_prior.npz")

NTHETA = 24          # 0..pi
NTAU = 36            # -pi..pi
SMOOTH = 1.0         # bins of Gaussian smoothing, wrapped in tau only
FLOOR = 1e-5


def pseudo_angles(ca: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """``(B, n-2)`` pseudo-bond-angles and ``(B, n-3)`` pseudo-torsions, radians."""
    x = np.atleast_3d(np.asarray(ca, float))
    if x.ndim == 2:
        x = x[None]
    v = x[:, 1:] - x[:, :-1]
    nv = np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-9)
    u = v / nv
    cos_t = np.clip(-(u[:, :-1] * u[:, 1:]).sum(-1), -1.0, 1.0)
    theta = np.arccos(cos_t)
    b1, b2, b3 = v[:, :-2], v[:, 1:-1], v[:, 2:]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m = np.cross(n1, b2 / np.maximum(np.linalg.norm(b2, axis=-1, keepdims=True), 1e-9))
    # Negated to match the IUPAC sign convention: a right-handed alpha helix must give
    # tau = +50 deg. Without this the raw formula on these difference vectors returns
    # -50 deg, which would make every handedness diagnostic read backwards. The learned
    # table is self-consistent either way, so this is a readability fix, not a scoring one
    # -- but a sign convention that disagrees with the literature is a trap.
    tau = -np.arctan2((m * n2).sum(-1), (n1 * n2).sum(-1))
    return theta, tau


def _bin(theta: np.ndarray, tau: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ti = np.clip((theta / np.pi * NTHETA).astype(int), 0, NTHETA - 1)
    ui = np.clip(((tau + np.pi) / (2 * np.pi) * NTAU).astype(int), 0, NTAU - 1)
    return ti, ui


def _smooth(h: np.ndarray) -> np.ndarray:
    """Gaussian smoothing: wrapped along tau, clamped along theta."""
    k = np.exp(-((np.arange(-2, 3)) ** 2) / (2 * SMOOTH ** 2))
    k = k / k.sum()
    out = np.zeros_like(h)
    for s, w in zip(range(-2, 3), k):
        out += w * np.roll(h, s, axis=1)
    h2 = np.zeros_like(out)
    for s, w in zip(range(-2, 3), k):
        h2 += w * np.take(out, np.clip(np.arange(NTHETA) - s, 0, NTHETA - 1), axis=0)
    return h2


def build(exclude_seq: str = "", entries: Optional[Sequence] = None) -> np.ndarray:
    """``(NTHETA, NTAU)`` negative log-probability table from held-out structures."""
    entries = entries if entries is not None else (
        pdb.holdout(exclude_seq) if exclude_seq else pdb.load())
    h = np.zeros((NTHETA, NTAU))
    for p in entries:
        if p.n < 4:
            continue
        th, tu = pseudo_angles(p.ca[None])
        ti, ui = _bin(th[0][:-1], tu[0])          # align theta_i with tau_i
        np.add.at(h, (ti, ui), 1.0)
    h = _smooth(h) + FLOOR * max(h.sum(), 1.0) / h.size
    return -np.log(h / h.sum())


@lru_cache(maxsize=8)
def table(exclude_seq: str = "") -> np.ndarray:
    if not exclude_seq and os.path.exists(CACHE):
        return np.load(CACHE)["nlp"]
    nlp = build(exclude_seq)
    if not exclude_seq:
        np.savez(CACHE, nlp=nlp)
    return nlp


class CATracePrior:
    """Batched ``-log P(theta, tau)`` per structure, and a mirror diagnostic."""

    def __init__(self, exclude_seq: str = ""):
        self.nlp = table(exclude_seq)

    def score(self, ca: np.ndarray) -> np.ndarray:
        """``(B,)`` mean negative log-probability over the trace. Lower is better."""
        arr = np.asarray(ca, float)
        single = arr.ndim == 2
        if single:
            arr = arr[None]
        if arr.shape[1] < 4:
            return np.zeros(len(arr))
        th, tu = pseudo_angles(arr)
        ti, ui = _bin(th[:, :-1], tu)
        out = self.nlp[ti, ui].mean(1)
        return float(out[0]) if single else out

    def handedness(self, ca: np.ndarray) -> np.ndarray:
        """``(B,)`` mean sign of the pseudo-torsion in helical windows.

        Positive is right-handed. A reflected structure returns the negation, which is what
        makes this a direct mirror detector rather than an indirect one.
        """
        arr = np.asarray(ca, float)
        if arr.ndim == 2:
            arr = arr[None]
        if arr.shape[1] < 4:
            return np.zeros(len(arr))
        th, tu = pseudo_angles(arr)
        helical = th[:, :-1] < np.radians(110.0)
        s = np.where(helical, np.sign(tu), 0.0)   # +1 right-handed, -1 left-handed
        n = np.maximum(helical.sum(1), 1)
        return s.sum(1) / n

    def mirror_penalty(self, ca: np.ndarray) -> np.ndarray:
        """``(B,)`` score difference between the structure and its mirror.

        Zero would mean the term cannot tell them apart; it is reported so the claim that
        this fixes the reflection blind spot is checkable rather than asserted.
        """
        arr = np.asarray(ca, float)
        if arr.ndim == 2:
            arr = arr[None]
        m = arr.copy()
        m[..., 2] *= -1.0
        return self.score(m) - self.score(arr)
