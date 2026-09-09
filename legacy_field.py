"""The Legacy (non-all-atom, knowledge-based) energy model as a BATCHED GENERATION field.

`hamiltonian.FoldingHamiltonian` evaluates the same physics one bitstring at a time, which
is fast enough to rescore a pool and far too slow to *drive* a search: `objective.
FoldObjective` scores 384 structures per VQE iteration, 60 iterations, 4 restarts, two
resolutions -- roughly 200k structures per fold. This module recomputes every coupled term
of `energy_terms` over a whole ``(B, n)`` state matrix with the same functional forms, so
the Legacy model can be blended into the objective the VQE actually minimises rather than
only used after the fact.

Where the forms differ from `energy_terms`, and why:

``hbond``     the scalar model runs a stable-sorted greedy one-donor-one-acceptor match.
              Here the same greedy rule is run as `n` vectorised rounds of "take the
              global minimum still available", which is the identical algorithm with a
              different tie-break. Ties are manufactured only by the `HB_E_FLOOR` clamp at
              -4.0, i.e. by interpenetrating pairs, which are rejected downstream anyway.
``aromatic``  the CB-CB fallback well only. The search representation is built with
              `chi_bits=False`, so no ring geometry exists on this path and the scalar
              model would take the same branch.
``steric``    identical, reusing `energy_terms._steric_layout` so the pair list, radii and
              N/O exemption cannot drift.

`verify` checks term-by-term agreement against `energy_terms.energy_components`.

Two things this module adds that the scalar model does not have, both of which are what
make it usable as a *generation* objective rather than a rescorer:

`LegacyField.standardise` fixes the term scales against a sample of random structures
drawn from the target's own state library -- never from the native, and never from any
other target -- so a weighted sum of eleven terms with different units and different
per-target dynamic ranges can be added to the distogram score at one transferable weight.

`FITTED_WEIGHTS` is a per-term weighting fitted to rank REAL candidate pools, not decoys.
See `work/legacy_generation.md`; `work/decoy_bank_invalid.md` explains why a synthetic
decoy bank cannot be used to fit a physics-based scorer at all.
"""
from typing import Dict, Optional, Sequence

import numpy as np

import energy_terms as et
import protein_geometry as geo

#: Term order. Same names as `energy_terms.TERM_NAMES` so weights are interchangeable.
TERMS = list(et.TERM_NAMES)

#: Named term families, for `LegacyField(use=...)`. `all` is the shipped Legacy model;
#: the rest isolate the families that could plausibly be independent of a CA-CA distance
#: prior. `hb` is pure backbone geometry with no amino-acid identity in it at all, which
#: makes it the cleanest test of independence: nothing in it can be a restatement of what
#: a sequence-conditioned distance predictor already knows.
SUBSETS = {
    "all": None,
    "hb": ("hbond_local", "hbond_longrange", "coop_helix", "coop_sheet"),
    "burial": ("contact", "solvation"),
    "packing": ("steric", "contact", "solvation", "compactness"),
}


class BatchLegacy:
    """Per-sequence tables built once; `terms(states)` scores a whole batch."""

    def __init__(self, sequence: str, representation, use_corrected_mj: bool = True):
        self.seq = sequence
        self.rep = representation
        self.n = n = len(sequence)
        self.burial, self.q, self.mj = et.sequence_arrays(sequence, use_corrected_mj)
        di, dj, sep = et.pair_index(n)
        self.di, self.dj, self.sep = di, dj, sep.astype(float)
        self.mj_pair = self.mj[di, dj]
        self.qq = self.q[di] * self.q[dj]
        self.elec_mask = (sep >= 2) & (self.qq != 0.0)
        self.m3 = sep >= 3
        self.rows = np.arange(n)
        self.idx_sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        # steric layout, shared with the scalar model
        self.st_names = ("N", "CA", "C", "O", "CB")
        self.st_ii, self.st_jj, self.st_lim = et._steric_layout(
            self.st_names, n, (), 2)
        # aromatic CB pairs at |i-j| >= 3
        idx = et.aromatic_indices(sequence)
        pairs = [(i, j) for a, i in enumerate(idx) for j in idx[a + 1:] if j - i >= 3]
        self.ar_i = np.array([p[0] for p in pairs], int)
        self.ar_j = np.array([p[1] for p in pairs], int)
        # rama table: per-residue, per-state, exact for a fixed library
        self.rama = np.array(
            [[et.rama_penalty(sequence[i], representation._phi[i, s],
                              representation._psi[i, s])
              for s in range(representation.n_states)] for i in range(n)])
        self.rg_target = 2.2 * (n ** 0.38)

    # ---------------------------------------------------------------- pieces
    def _steric(self, coords) -> np.ndarray:
        A = np.concatenate([coords[k] for k in self.st_names], axis=1)   # (B, 5n, 3)
        d = np.linalg.norm(A[:, self.st_ii] - A[:, self.st_jj], axis=2)
        over = np.maximum(0.0, self.st_lim[None, :] - d)
        return (over ** 2).sum(1)

    def _hbonds(self, coords):
        """Greedy one-donor-one-acceptor DSSP match, batched.

        Returns ``(local, longrange, bond)`` with `bond` the (B, n, n) boolean matrix of
        matched (donor, acceptor) pairs -- what the cooperativity terms read.
        """
        N, C, O = coords["N"], coords["C"], coords["O"]
        B, n = N.shape[0], self.n
        H = np.full((B, n, 3), np.nan)
        if n >= 2:
            d = C[:, :-1] - O[:, :-1]
            nd = np.linalg.norm(d, axis=2)
            ok0 = nd > 1e-6
            sh = N[:, 1:] + np.divide(d, nd[..., None], out=np.zeros_like(d),
                                      where=ok0[..., None])
            H[:, 1:][ok0] = sh[ok0]
        valid = np.isfinite(H).all(axis=2)
        dON = np.linalg.norm(N[:, :, None, :] - O[:, None, :, :], axis=3)
        dCH = np.linalg.norm(C[:, None, :, :] - H[:, :, None, :], axis=3)
        dOH = np.linalg.norm(H[:, :, None, :] - O[:, None, :, :], axis=3)
        dCN = np.linalg.norm(N[:, :, None, :] - C[:, None, :, :], axis=3)
        with np.errstate(divide="ignore", invalid="ignore"):
            E = 0.084 * 332.0 * (1.0 / dON + 1.0 / dCH - 1.0 / dOH - 1.0 / dCN)
        E = np.maximum(np.nan_to_num(E, nan=0.0, posinf=0.0, neginf=geo.HB_E_FLOOR),
                       geo.HB_E_FLOOR)
        ok = valid[:, :, None] & (self.idx_sep[None] >= 2)
        ok &= (dON > geo.HB_MIN_ON) & (dCH > 0.5) & (dOH > 0.5) & (dCN > 0.5)
        ok &= np.isfinite(E) & (E < -0.5)

        # greedy: repeatedly take the lowest still-admissible pair, block its donor and
        # acceptor. `n` rounds suffice -- each round consumes one donor.
        big = 1e9
        avail = ok.copy()
        bond = np.zeros_like(ok)
        Ef = np.where(ok, E, big)
        bi = np.arange(B)
        for _ in range(n):
            flat = np.where(avail, Ef, big).reshape(B, n * n)
            k = flat.argmin(1)
            v = flat[bi, k]
            live = v < big
            if not live.any():
                break
            d_i, a_j = k // n, k % n
            bond[bi[live], d_i[live], a_j[live]] = True
            avail[bi[live], d_i[live], :] = False
            avail[bi[live], :, a_j[live]] = False
        # desolvation cost of +1.0, matching `energy_terms.hbond_terms`
        e = np.where(bond, E + 1.0, 0.0)
        e = np.where(e < 0.0, e, 0.0)
        bond = bond & (E + 1.0 < 0.0)
        sepm = self.idx_sep[None]
        local = np.where(sepm < et.HB_LONGRANGE_SEP, e, 0.0).sum((1, 2))
        lr = np.where(sepm >= et.HB_LONGRANGE_SEP, e, 0.0).sum((1, 2))
        return local, lr, bond

    def _coop(self, bond):
        """Consecutive n-turns and consecutive antiparallel ladder rungs, batched."""
        n = self.n
        d, a = np.indices((n, n))
        hel = bond & (d - a >= 3)[None] & (d - a <= 4)[None]
        run = np.zeros(bond.shape[0])
        if n >= 2:
            nxt = np.zeros_like(hel)
            nxt[:, :-1, :-1] = hel[:, 1:, 1:]
            run = (hel & nxt).sum((1, 2)).astype(float)
        helix = -np.minimum(run, et.COOP_RUN_CAP)
        lad = bond & (np.abs(d - a) >= 2)[None]
        nxt = np.zeros_like(lad)
        if n >= 3:
            nxt[:, :-2, 2:] = lad[:, 2:, :-2]
        sheet = -(lad & nxt).sum((1, 2)).astype(float)
        return helix, sheet

    # ---------------------------------------------------------------- public
    def terms_from_coords(self, coords: Dict[str, np.ndarray],
                          states: Optional[np.ndarray] = None,
                          phi: Optional[np.ndarray] = None,
                          psi: Optional[np.ndarray] = None) -> np.ndarray:
        """``(B, 11)`` unweighted term values, columns ordered as `TERMS`."""
        CA = coords["CA"]
        CB = coords.get("CB", CA)
        B, n = CA.shape[0], self.n
        dcb = np.linalg.norm(CB[:, self.di, :] - CB[:, self.dj, :], axis=2)
        out = np.zeros((B, len(TERMS)))
        col = {t: k for k, t in enumerate(TERMS)}

        out[:, col["steric"]] = self._steric(coords)
        out[:, col["contact"]] = (et.switch(dcb[:, self.m3], 4.5, 8.5)
                                  * self.mj_pair[None, self.m3]).sum(1)
        local, lr, bond = self._hbonds(coords)
        out[:, col["hbond_local"]] = local
        out[:, col["hbond_longrange"]] = lr
        h, s = self._coop(bond)
        out[:, col["coop_helix"]] = h
        out[:, col["coop_sheet"]] = s
        # solvation: burial-weighted coordination over CB pairs
        sw = et.switch(dcb, 6.0, 10.0)
        coord = np.zeros((B, n))
        np.add.at(coord.T, self.di, sw.T)
        np.add.at(coord.T, self.dj, sw.T)
        out[:, col["solvation"]] = -(coord @ (self.burial / et.BURIAL_NORM))
        # electrostatic
        if self.elec_mask.any():
            de = np.maximum(dcb[:, self.elec_mask], 2.0)
            dr = dcb[:, self.elec_mask]
            out[:, col["electrostatic"]] = (
                (et.COULOMB / et.DIELECTRIC)
                * (self.qq[None, self.elec_mask] / de) * np.exp(-dr / 8.0)).sum(1)
        # aromatic, CB proxy
        if len(self.ar_i):
            da = np.linalg.norm(CB[:, self.ar_i, :] - CB[:, self.ar_j, :], axis=2)
            out[:, col["aromatic"]] = -np.exp(
                -((da - et.AROM_CB_DIST) / et.AROM_CB_WIDTH) ** 2).sum(1)
        # torsion
        if states is not None:
            S = np.asarray(states, int)
            out[:, col["torsion"]] = self.rama[self.rows[None, :], S].sum(1)
        elif phi is not None:
            out[:, col["torsion"]] = np.array(
                [et.torsion_term(self.seq, phi[b], psi[b]) for b in range(B)])
        # compactness
        rg = np.sqrt(((CA - CA.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        out[:, col["compactness"]] = np.maximum(0.0, rg - self.rg_target) ** 2
        return out

    def terms(self, states: np.ndarray) -> np.ndarray:
        S = np.asarray(states, int)
        phi = self.rep._phi[self.rows[None, :], S]
        psi = self.rep._psi[self.rows[None, :], S]
        coords = geo.build_backbone_batch(phi, psi)
        return self.terms_from_coords(coords, states=S)

# Weighting
#: Per-term weights for the *generation* field. Fitted on real candidate pools of the
#: 24 development peptides (`work/pools`, dev split), never on a benchmark target and
#: never on the synthetic decoy bank. See `work/fit_legacy.py`.
FITTED_WEIGHTS: Dict[str, float] = dict(et.DEFAULT_WEIGHTS)


class LegacyField:
    """A single scalar per structure: standardised, weighted Legacy energy.

    `standardise` draws `n_ref` random state assignments from the target's own library and
    records each term's mean and spread on that population. The field returned is then
    dimensionless and comparable across targets and across chain lengths, which is what
    lets one blend weight `w_legacy` transfer. Nothing in the reference sample depends on
    the native structure.
    """

    def __init__(self, sequence: str, representation,
                 weights: Optional[Dict[str, float]] = None,
                 n_ref: int = 512, seed: int = 0,
                 use: Optional[Sequence[str]] = None, dist=None):
        self.bl = BatchLegacy(sequence, representation)
        w = dict(FITTED_WEIGHTS if weights is None else weights)
        self.w = np.array([w.get(t, 0.0) for t in TERMS], float)
        if use is not None:
            keep = set(use)
            self.w = np.array([self.w[k] if t in keep else 0.0
                               for k, t in enumerate(TERMS)])
        self.mu = np.zeros(len(TERMS))
        self.sd = np.ones(len(TERMS))
        self._standardise(representation, n_ref, seed, dist)

    def _standardise(self, rep, n_ref: int, seed: int, dist) -> None:
        rng = np.random.default_rng(seed)
        S = rng.integers(0, rep.n_states, size=(n_ref, rep.n_residues))
        rows = np.arange(rep.n_residues)
        coords = geo.build_backbone_batch(rep._phi[rows[None], S],
                                          rep._psi[rows[None], S])
        T = self.bl.terms_from_coords(coords, states=S)
        self.mu = T.mean(0)
        self.sd = np.maximum(T.std(0), 1e-6)
        # the weighted sum's own scale on the same population, so `w_legacy` is a
        # multiple of "one standard deviation of the random-structure population"
        z = (T - self.mu) / self.sd
        v = z @ self.w
        self.centre = float(v.mean())
        self.scale = float(max(v.std(), 1e-6))
        # ... and, when the distance term is supplied, expressed in ITS units on the same
        # sample. Without this `w_legacy` would mean something different on every target,
        # because the Bayes-risk score's spread depends on how confident the prior is.
        self.unit = 1.0
        if dist is not None:
            d = np.asarray(dist.score(coords["CA"]), float)
            self.unit = float(max(d.std(), 1e-9))

    def _reduce(self, T: np.ndarray) -> np.ndarray:
        z = (T - self.mu[None, :]) / self.sd[None, :]
        return self.unit * (z @ self.w - self.centre) / self.scale

    def score(self, states: np.ndarray) -> np.ndarray:
        return self._reduce(self.bl.terms(states))

    def score_from_coords(self, coords, states=None, phi=None, psi=None) -> np.ndarray:
        return self._reduce(self.bl.terms_from_coords(coords, states=states,
                                                      phi=phi, psi=psi))


def verify(sequence: str, representation, batch: int = 24, seed: int = 0):
    """Max per-term deviation of `BatchLegacy` from `energy_terms.energy_components`."""
    bl = BatchLegacy(sequence, representation)
    rng = np.random.default_rng(seed)
    S = rng.integers(0, representation.n_states, size=(batch, len(sequence)))
    T = bl.terms(S)
    rows = np.arange(len(sequence))
    worst = {t: 0.0 for t in TERMS}
    for k in range(batch):
        phi = representation._phi[rows, S[k]]
        psi = representation._psi[rows, S[k]]
        coords = geo.build_backbone(phi, psi)
        ref = et.energy_components(sequence, coords, phi, psi)
        for c, t in enumerate(TERMS):
            worst[t] = max(worst[t], abs(ref[t] - T[k, c]))
    return worst
