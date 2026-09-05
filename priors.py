"""Sequence-conditioned learned structural priors.

Two priors, both fitted on `peptide_db` with the target's sequence and its homologs held
out, and both consumed as *tables* so that scoring a structure during search is a gather
rather than a model call.

**Distance prior (distogram).** For every residue pair (i, j) of the target, a
distribution over the CA-CA distance, predicted from sequence alone. This is the piece the
old objective did not have. Every descriptor in `batch_features` is a global shape
statistic -- radius of gyration, contact counts, contact order -- and global statistics
cannot separate a 1.5 A structure from a 3 A one, because both have the same size and
roughly the same number of contacts; they differ in *which* residues are in contact. The
distogram scores exactly that.

Scoring is ``mean_ij -log P(bin(d_ij) | i, j, sequence)``, the negative log-likelihood of
the structure under the prior. The (n, n, NBINS) log-probability table is built once per
target, so a whole batch of structures is scored by digitising a distance matrix and
gathering -- microseconds, not model calls.

**Torsion-state prior (MRF).** Per-residue state probabilities and nearest-neighbour state
couplings over the discrete torsion library, estimated from the same held-out database by
projecting matching database residues onto the target's own per-residue library. Gives
three things: a local energy term, a warm start for the VQE (the ansatz's RY angles set the
per-qubit marginals directly), and a proposal distribution for refinement.

Neither prior ever sees the target structure. `peptide_db.holdout` removes the target
sequence and everything above 0.6 identity to it before a single parameter is fitted.
"""
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import peptide_db as pdb
import representations as reps

BASE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------- residue features
#: Chou-Fasman helix/sheet propensity, Kyte-Doolittle hydropathy, formal charge at pH 7,
#: and side-chain volume (A^3). Five numbers per residue type -- enough to let the model
#: generalise across amino acids instead of memorising 20 identities from 780 peptides.
_PROPS = {
    #      helix sheet  hydro  charge  volume
    "A": (1.42, 0.83,  1.8,  0.0,  88.6), "C": (0.70, 1.19,  2.5,  0.0, 108.5),
    "D": (1.01, 0.54, -3.5, -1.0, 111.1), "E": (1.51, 0.37, -3.5, -1.0, 138.4),
    "F": (1.13, 1.38,  2.8,  0.0, 189.9), "G": (0.57, 0.75, -0.4,  0.0,  60.1),
    "H": (1.00, 0.87, -3.2,  0.1, 153.2), "I": (1.08, 1.60,  4.5,  0.0, 166.7),
    "K": (1.16, 0.74, -3.9,  1.0, 168.6), "L": (1.21, 1.30,  3.8,  0.0, 166.7),
    "M": (1.45, 1.05,  1.9,  0.0, 162.9), "N": (0.67, 0.89, -3.5,  0.0, 114.1),
    "P": (0.57, 0.55, -1.6,  0.0, 112.7), "Q": (1.11, 1.10, -3.5,  0.0, 143.8),
    "R": (0.98, 0.93, -4.5,  1.0, 173.4), "S": (0.77, 0.75, -0.8,  0.0,  89.0),
    "T": (0.83, 1.19, -0.7,  0.0, 116.1), "V": (1.06, 1.70,  4.2,  0.0, 140.0),
    "W": (1.08, 1.37, -0.9,  0.0, 227.8), "Y": (0.69, 1.47, -1.3,  0.0, 193.6),
}
_DEFAULT = tuple(float(np.mean([v[k] for v in _PROPS.values()])) for k in range(5))


def _prop(seq: str) -> np.ndarray:
    return np.array([_PROPS.get(a, _DEFAULT) for a in seq], float)


#: Distance bin edges, A. Fine where CA-CA distances are structurally informative
#: (5-14 A spans an i,i+3 helical turn through a hairpin cross-strand pair) and coarse
#: beyond, where "far apart" is all the information there is.
BIN_EDGES = np.array([4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0,
                      11.0, 12.5, 14.0, 16.0, 19.0, 23.0])
NBINS = len(BIN_EDGES) + 1


def pair_features(seq: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(features, i, j)`` for every pair with ``j - i >= 2`` of one sequence."""
    n = len(seq)
    P = _prop(seq)
    # Local windows: mean propensity over +-2 residues, the classic secondary-structure
    # signal, computed by stacking shifted copies so it costs nothing.
    pad = np.vstack([P[:1], P[:1], P, P[-1:], P[-1:]])
    W = np.stack([pad[k:k + n] for k in range(5)]).mean(0)
    i, j = np.triu_indices(n, k=2)
    s = (j - i).astype(float)
    between = np.stack([P[a + 1:b].mean(0) if b > a + 1 else P[a]
                        for a, b in zip(i, j)])
    f = np.column_stack([
        s, np.log(s), s / n, np.full(i.shape, float(n)),
        i / n, j / n, np.minimum(i, n - 1 - j).astype(float),
        P[i], P[j], W[i], W[j],
        P[i] * P[j], np.abs(P[i] - P[j]),
        between,
    ])
    return f.astype(np.float32), i, j


def fit_distogram_model(entries: Sequence[pdb.Peptide], seed: int = 0,
                        max_iter: int = 200):
    from sklearn.ensemble import HistGradientBoostingClassifier
    X, Y = [], []
    for p in entries:
        f, i, j = pair_features(p.seq)
        if not len(f):
            continue
        d = np.linalg.norm(p.ca[i] - p.ca[j], axis=1)
        X.append(f)
        Y.append(np.digitize(d, BIN_EDGES))
    X = np.vstack(X)
    Y = np.concatenate(Y)
    m = HistGradientBoostingClassifier(max_iter=max_iter, learning_rate=0.1, max_depth=7,
                                       l2_regularization=1.0, random_state=seed,
                                       early_stopping=False)
    m.fit(X, Y)
    return m


class DistogramPrior:
    """Per-target log P(distance bin) table, shape ``(n, n, NBINS)``."""

    def __init__(self, sequence: str, logp: np.ndarray):
        self.seq = sequence
        self.n = len(sequence)
        self.logp = logp
        self._iu = np.triu_indices(self.n, k=2)
        self._flat = logp[self._iu[0], self._iu[1]]          # (npairs, NBINS)

    @classmethod
    def fit(cls, sequence: str, entries: Optional[Sequence[pdb.Peptide]] = None,
            model=None, seed: int = 0) -> "DistogramPrior":
        if model is None:
            model = fit_distogram_model(
                entries if entries is not None else pdb.holdout(sequence), seed)
        n = len(sequence)
        f, i, j = pair_features(sequence)
        prob = model.predict_proba(f)
        cols = np.asarray(model.classes_, int)
        tab = np.full((len(f), NBINS), 1e-4)
        tab[:, cols] = np.maximum(prob, 1e-4)
        tab /= tab.sum(1, keepdims=True)
        full = np.full((n, n, NBINS), 1.0 / NBINS)
        full[i, j] = tab
        full[j, i] = tab
        return cls(sequence, np.log(full).astype(np.float32))

    def score(self, ca: np.ndarray) -> np.ndarray:
        """Mean per-pair negative log-likelihood for a batch ``(B, n, 3)`` of CA traces."""
        arr = np.asarray(ca, float)
        single = arr.ndim == 2
        if single:
            arr = arr[None]
        i, j = self._iu
        d = np.linalg.norm(arr[:, i, :] - arr[:, j, :], axis=-1)
        b = np.searchsorted(BIN_EDGES, d)
        lp = np.take_along_axis(self._flat[None], b[:, :, None], axis=2)[:, :, 0]
        out = -lp.mean(1)
        return float(out[0]) if single else out

    def expected_distances(self) -> np.ndarray:
        """Posterior-mean CA-CA distance per pair, for diagnostics and restraints."""
        centres = np.concatenate([[BIN_EDGES[0] - 0.5],
                                  0.5 * (BIN_EDGES[1:] + BIN_EDGES[:-1]),
                                  [BIN_EDGES[-1] + 2.0]])
        return (np.exp(self.logp) * centres).sum(-1)


# --------------------------------------------------------------------- torsion MRF
class TorsionMRF:
    """Per-residue state marginals plus nearest-neighbour state couplings.

    Estimated by projecting held-out database residues onto the target's per-residue
    library: for slot ``i``, every database residue of the same type at a comparable
    fractional chain position contributes its nearest state. Couplings use consecutive
    database pairs whose residue types match ``(seq[i], seq[i+1])``.
    """

    def __init__(self, h: np.ndarray, J: np.ndarray):
        self.h = h                # (n, k)      -log marginal
        self.J = J                # (n-1, k, k) -log (joint / marginal product)

    @classmethod
    def fit(cls, sequence: str, table: np.ndarray,
            entries: Sequence[pdb.Peptide], smooth: float = 1.0,
            pos_tol: float = 0.34) -> "TorsionMRF":
        n, k, _ = table.shape
        counts = np.full((n, k), smooth)
        pair = np.full((n - 1, k, k), smooth / k)
        cphi, sphi = np.cos(table[:, :, 0]), np.sin(table[:, :, 0])
        cpsi, spsi = np.cos(table[:, :, 1]), np.sin(table[:, :, 1])

        def nearest(i: int, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
            d = ((np.cos(phi)[:, None] - cphi[i]) ** 2
                 + (np.sin(phi)[:, None] - sphi[i]) ** 2
                 + (np.cos(psi)[:, None] - cpsi[i]) ** 2
                 + (np.sin(psi)[:, None] - spsi[i]) ** 2)
            return d.argmin(1)

        by_aa: Dict[str, List[Tuple[float, float, float, str, float, float]]] = {}
        for p in entries:
            for t in range(p.n):
                by_aa.setdefault(p.seq[t], []).append(
                    (t / max(1, p.n - 1), float(p.phi[t]), float(p.psi[t]),
                     p.seq[t + 1] if t + 1 < p.n else "",
                     float(p.phi[t + 1]) if t + 1 < p.n else np.nan,
                     float(p.psi[t + 1]) if t + 1 < p.n else np.nan))

        for i in range(n):
            frac = i / max(1, n - 1)
            rows = [r for r in by_aa.get(sequence[i], ())
                    if abs(r[0] - frac) <= pos_tol]
            if rows:
                s = nearest(i, np.array([r[1] for r in rows]),
                            np.array([r[2] for r in rows]))
                np.add.at(counts[i], s, 1.0)
            if i < n - 1:
                nxt = [r for r in rows if r[3] == sequence[i + 1]]
                if nxt:
                    s0 = nearest(i, np.array([r[1] for r in nxt]),
                                 np.array([r[2] for r in nxt]))
                    s1 = nearest(i + 1, np.array([r[4] for r in nxt]),
                                 np.array([r[5] for r in nxt]))
                    np.add.at(pair[i], (s0, s1), 1.0)

        marg = counts / counts.sum(1, keepdims=True)
        h = -np.log(marg)
        J = np.zeros((max(0, n - 1), k, k))
        for i in range(n - 1):
            pj = pair[i] / pair[i].sum()
            J[i] = (-np.log(np.maximum(pj, 1e-12))
                    + np.log(marg[i])[:, None] + np.log(marg[i + 1])[None, :])
        return cls(h.astype(np.float32), J.astype(np.float32))

    def score(self, states: np.ndarray) -> np.ndarray:
        """``states`` (B, n) of state indices -> (B,) mean per-residue MRF energy."""
        S = np.atleast_2d(np.asarray(states, int))
        n = S.shape[1]
        e = self.h[np.arange(n)[None], S].sum(1)
        if n > 1 and len(self.J):
            e = e + self.J[np.arange(n - 1)[None], S[:, :-1], S[:, 1:]].sum(1)
        return e / n

    def marginals(self) -> np.ndarray:
        p = np.exp(-self.h)
        return p / p.sum(1, keepdims=True)
