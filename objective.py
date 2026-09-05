"""The composite folding objective the VQE minimises, evaluated in batch.

Three terms, each earning its place from a measurement rather than from physical
plausibility:

``distogram``  weighted L1 Bayes risk against the predicted CA-CA distance distribution.
               This is the only term with real near-native discrimination; see
               `distogram` for the measurements that put it at the centre.
``mrf``        the local torsion prior: per-residue state marginals plus nearest-neighbour
               couplings, estimated from the held-out database. It is a *local* term, so
               on its own it has no fold in it -- its job is to keep the chain in
               populated Ramachandran basins and to supply the VQE's warm start, not to
               choose a tertiary arrangement.
``clash``      excluded volume on CA and CB. A feasibility filter, not a discriminator:
               flat across everything the search should be choosing between, and steeply
               positive on the structures it must not return.
``legacy``     the coarse knowledge-based potential of `energy_terms`, batched by
               `legacy_field`, standardised against a random sample from this target's own
               state library and blended at `w_legacy` *in units of the distance term's own
               spread on that same sample*, so one weight transfers across targets and
               chain lengths. Off by default -- `work/legacy_generation.md` has the
               measurement that decides the weight, and it is not a large number.

Everything is a function of the ``(B, n)`` state matrix, and every step from states to
score is batched, so one call scores a whole VQE shot budget.
"""
from typing import Dict, Optional

import numpy as np

import protein_geometry as geo


class FoldObjective:
    """Batched score over discrete torsion-state assignments."""

    def __init__(self, sequence: str, rep, dist=None, mrf=None, catr=None,
                 w_dist: float = 1.0, w_mrf: float = 0.35, w_clash: float = 1.0,
                 w_catrace: float = 0.0, legacy=None, w_legacy: float = 0.0,
                 clash_ca: float = 4.1, clash_cb: float = 3.6):
        self._cons_d = None
        self._cons_w = None
        self.w_cons = 0.0
        self.seq = sequence
        self.rep = rep
        self.n = rep.n_residues
        self.dist = dist
        self.mrf = mrf
        # CA-trace (theta, tau) prior. The distance term is invariant under reflection, so
        # without this the objective literally cannot tell a structure from its mirror --
        # and 15-24% of the pool on 1LE0 is mirrored, with the mirror selected on 2 of 3
        # seeds. See `catrace`.
        self.catr = catr
        self.w_catrace = w_catrace
        # Coarse physical potential as a GENERATION term. `legacy_field.LegacyField`
        # standardises itself on a random sample from `rep`, so `w_legacy` is a multiple of
        # one standard deviation of that population -- and `match_scale` puts it in the
        # distance term's units, so the same number means the same thing on every target.
        self.legacy = legacy
        self.w_legacy = float(w_legacy)
        self.w_dist, self.w_mrf, self.w_clash = w_dist, w_mrf, w_clash
        self.clash_ca, self.clash_cb = clash_ca, clash_cb
        self.rows = np.arange(self.n)
        self._iu = np.triu_indices(self.n, k=3)
        self.n_calls = 0
        self.n_structures = 0

    # -- geometry ---------------------------------------------------------
    def build(self, states: np.ndarray) -> Dict[str, np.ndarray]:
        S = np.asarray(states, int)
        return geo.build_backbone_batch(self.rep._phi[self.rows[None], S],
                                        self.rep._psi[self.rows[None], S])

    def clash(self, coords: Dict[str, np.ndarray]) -> np.ndarray:
        i, j = self._iu
        out = np.zeros(len(coords["CA"]))
        for key, cut in (("CA", self.clash_ca), ("CB", self.clash_cb)):
            if key not in coords:
                continue
            d = np.linalg.norm(coords[key][:, i] - coords[key][:, j], axis=-1)
            out += np.maximum(cut - d, 0.0).sum(1)
        return out

    def set_consensus(self, distances: Optional[np.ndarray],
                      weights: Optional[np.ndarray] = None,
                      w_cons: float = 0.0) -> None:
        """Add an L1 restraint toward a consensus distance matrix from a first pass.

        The search's own candidate pool carries distance information the prior does not:
        the pool's best member is far closer to native than the one the prior selects. A
        second pass restrained toward the pool's score-weighted mean distances is therefore
        not merely echoing the prior. It is also the obvious way to amplify the prior's
        mistakes, so `w_cons` is small and is fitted on the development peptides.

        `distances` is indexed like `dist.i` / `dist.j`; pass None to clear.
        """
        self._cons_d = None if distances is None else np.asarray(distances, float)
        self._cons_w = (None if weights is None else
                        np.asarray(weights, float) / max(np.mean(weights), 1e-12))
        self.w_cons = float(w_cons)

    def _consensus(self, ca: np.ndarray) -> np.ndarray:
        d = self.dist
        dd = np.linalg.norm(ca[:, d.i, :] - ca[:, d.j, :], axis=-1)
        err = np.abs(dd - self._cons_d[None, :])
        if self._cons_w is not None:
            err = err * self._cons_w[None, :]
        return err.mean(1)

    # -- score ------------------------------------------------------------
    def score_from(self, states: np.ndarray, coords: Dict[str, np.ndarray]) -> np.ndarray:
        e = np.zeros(len(coords["CA"]))
        if self.dist is not None and self.w_dist:
            e += self.w_dist * self.dist.score(coords["CA"])
        if self.mrf is not None and self.w_mrf:
            e += self.w_mrf * self.mrf.score(states)
        if self.w_clash:
            e += self.w_clash * self.clash(coords)
        if self.catr is not None and self.w_catrace:
            e += self.w_catrace * self.catr.score(coords["CA"])
        if self.legacy is not None and self.w_legacy:
            e += self.w_legacy * self.legacy.score_from_coords(coords, states=states)
        if self.w_cons and self._cons_d is not None:
            e += self.w_cons * self._consensus(coords["CA"])
        return e

    def __call__(self, states: np.ndarray):
        """``states`` (B, n) -> ``(scores (B,), CA (B, n, 3))``."""
        coords = self.build(states)
        self.n_calls += 1
        self.n_structures += len(coords["CA"])
        return self.score_from(states, coords), coords["CA"]

    # -- continuous form, for refinement ----------------------------------
    def score_angles(self, phi: np.ndarray, psi: np.ndarray):
        coords = geo.build_backbone_batch(phi, psi)
        e = np.zeros(len(coords["CA"]))
        if self.dist is not None and self.w_dist:
            e += self.w_dist * self.dist.score(coords["CA"])
        if self.w_clash:
            e += self.w_clash * self.clash(coords)
        if self.catr is not None and self.w_catrace:
            e += self.w_catrace * self.catr.score(coords["CA"])
        if self.legacy is not None and self.w_legacy:
            e += self.w_legacy * self.legacy.score_from_coords(coords, phi=phi, psi=psi)
        if self.w_cons and self._cons_d is not None:
            e += self.w_cons * self._consensus(coords["CA"])
        self.n_calls += 1
        self.n_structures += len(coords["CA"])
        return e, coords["CA"]
