"""Batched continuous Legacy energy, and local continuous refinement under it.

Two things live here, and they exist for the same reason.

`BatchedLegacy` evaluates `energy_terms.energy_components` for a whole batch of (phi, psi)
arrays at once. Profiled on a 10-mer (`work/prof_legacy.py`), the scalar entry point costs
~320 us per structure, of which the DSSP H-bond matrix is ~90 us, the soft-sphere steric
term ~28 us, and the Ramachandran sum a further ~40 us once the torsions are continuous
and its `lru_cache` stops hitting. That is fine for scoring a few thousand bitstrings; it
is not fine for a local optimiser needing O(4n) evaluations per finite-difference
gradient. Batched it is 3.4-9.1x faster (~5x typical, `work/speed2.py`), and everything
here is numerically equivalent to the scalar path -- `selftest` holds every term to a
relative 1e-9 and measures ~1e-13 -- differing only in floating-point summation order.

`refine` is the reason the batched form was needed: a bounded local relaxation of the
continuous torsions under the Legacy energy, starting from a discrete encoded structure.
The Legacy objective is NOT smooth -- the H-bond term admits a pair only below a hard
energy cutoff, and the two cooperativity terms are integer-valued run counts -- so a pure
quasi-Newton descent on a finite-difference gradient stalls on the plateaus between
integer changes. The optimiser is therefore a deterministic batched pattern search
(isotropic proposals with a shrinking radius, hard-clipped to a box around the start)
with an optional L-BFGS-B polish on the smooth part. Both are deterministic given the
seed.

Nothing here mutates `energy_terms`; `DEFAULT_WEIGHTS` behaviour is reproduced exactly
when `weights` is None.

WHAT THE REFINER WAS MEASURED TO DO, so nobody re-runs this hoping for a different answer
(`work/refine_study.py`, `work/refine_native.py`, 12 benchmark targets, k=8):

    start                       box 5 deg   box 15 deg   box 30 deg
    oracle floor states         +0.27 A      +0.65 A      +1.43 A     (0/12 improved
    native torsions             +0.42 A      +0.90 A      +1.84 A      at 15 and 30)
    lowest-Legacy-energy decoy  -0.01 A      +0.12 A      +0.05 A

The energy falls hard in every one of those cells (native: -6 to -11 kcal), so this is not
an optimiser failure -- the Legacy minima are simply 1-2 A away from the native basin, and
the native is not even a local minimum of the model. The decoy-start row is neutral
because it starts ~3.9 A out where CA-RMSD is flat, not because refinement is working.
Reweighting does not rescue it: the weights fitted in `work/reweight.py` (which do improve
in-band RANKING out of sample) leave the drift unchanged or slightly worse. Use this
module to measure the objective, not to improve a structure.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import energy_terms as et
import protein_geometry as geo

__all__ = ["BatchedLegacy", "batched_energy", "refine", "refine_angles",
           "weights_from_vector", "vector_from_weights", "selftest",
           "DEFAULT_BOX_DEG", "DEFAULT_STEPS"]

DEFAULT_BOX_DEG = 15.0
DEFAULT_STEPS = 60

_ATOM_ORDER = ("N", "CA", "C", "O", "CB")
_HB_PREFACTOR = 0.084 * 332.0


# ==========================================================================
# Batched energy
# ==========================================================================
class BatchedLegacy:
    """Legacy energy for a batch of structures of one fixed sequence.

    Every layout that depends only on the sequence -- pair lists, MJ values, charge
    products, steric radii, the Ramachandran basin weight table -- is built once in the
    constructor, so a call costs geometry plus arithmetic and no Python per-residue work
    except the H-bond greedy match (which sees ~1 admissible pair per structure).

    The aromatic term uses the CB fallback, which is what `energy_terms.aromatic_term`
    itself does whenever `rings` is None -- i.e. for every representation without chi1
    bits, which is the case for the continuous refinement this module exists to serve.
    Pass `rings` to the scalar path if you need real ring geometry; this class does not
    model it and `selftest` compares against the rings-free scalar call.
    """

    def __init__(self, sequence: str, weights: Optional[Dict[str, float]] = None,
                 use_corrected_mj: bool = True):
        self.sequence = sequence.strip().upper()
        self.n = len(self.sequence)
        self.weights = dict(et.DEFAULT_WEIGHTS if weights is None else weights)
        self.use_corrected_mj = bool(use_corrected_mj)

        n = self.n
        burial, q, mj = et.sequence_arrays(self.sequence, self.use_corrected_mj)
        di, dj, sep = et.pair_index(n)
        self._di, self._dj, self._sep = di, dj, sep

        # -- contact: MJ value per pair, restricted to |i-j| >= 3 ------------
        self._m3 = sep >= 3
        self._mj3 = mj[di[self._m3], dj[self._m3]]

        # -- solvation: each pair contributes its switch to both endpoints ---
        self._sol_w = (burial[di] + burial[dj]) / et.BURIAL_NORM

        # -- electrostatic: only charged pairs at |i-j| >= 2 -----------------
        qq = q[di] * q[dj]
        self._el_mask = (sep >= 2) & (qq != 0.0)
        self._el_qq = qq[self._el_mask]

        # -- aromatic (CB fallback) pairs ------------------------------------
        idx = et.aromatic_indices(self.sequence)
        pairs = [(i, j) for a, i in enumerate(idx) for j in idx[a + 1:]
                 if j - i >= 3]
        self._arom_i = np.array([p[0] for p in pairs], dtype=int)
        self._arom_j = np.array([p[1] for p in pairs], dtype=int)

        # -- steric layout (no rings) ----------------------------------------
        # The scalar term gathers two (P, 3) coordinate blocks per structure and takes a
        # norm; batched, that gather is the single most expensive operation in the model
        # (P ~ 830 pairs even for a 10-mer). Instead the limits are scattered into a dense
        # (A, A) matrix, so the whole term is one batched Gram product plus a comparison
        # on the squared distances. Pairs that are not in the list get limit 0, and no
        # squared distance is below 0, so they can never contribute.
        self._st_ii, self._st_jj, self._st_lim = et._steric_layout(
            _ATOM_ORDER, n, (), 2)
        A = len(_ATOM_ORDER) * n
        self._st_A = A
        self._st_lim_mat = np.zeros((A, A))
        self._st_lim_mat[self._st_ii, self._st_jj] = self._st_lim

        # -- Ramachandran basin weights --------------------------------------
        self._rama_W = np.zeros((n, len(et._RAMA_BASINS)))
        for i, aa in enumerate(self.sequence):
            for k, (_, _, _, depth) in enumerate(et._RAMA_BASINS):
                w = depth
                if k == 0 and aa in et._HELIX_FORMERS:
                    w += 0.5
                if k == 1 and aa in et._SHEET_FORMERS:
                    w += 0.5
                if k == 3 and aa != "G":
                    w *= 0.2
                self._rama_W[i, k] = w
        self._rama_pc = np.array([b[0] for b in et._RAMA_BASINS])
        self._rama_sc = np.array([b[1] for b in et._RAMA_BASINS])
        self._rama_two_sig2 = np.array([2.0 * b[2] ** 2 for b in et._RAMA_BASINS])
        self._is_pro = np.array([a == "P" for a in self.sequence])
        self._gly_shift = -0.2 * float(sum(1 for a in self.sequence if a == "G"))

        # -- H-bond bookkeeping -----------------------------------------------
        ii = np.arange(n)
        self._hb_sep_ok = np.abs(ii[:, None] - ii[None, :]) >= 2

        self._rg_target = 2.2 * (n ** 0.38)

    # -- geometry ----------------------------------------------------------
    def build(self, phi: np.ndarray, psi: np.ndarray) -> Dict[str, np.ndarray]:
        return geo.build_backbone_batch(np.atleast_2d(phi), np.atleast_2d(psi))

    # -- individual terms ---------------------------------------------------
    @staticmethod
    def _sqdist(A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """``(B, na, nb)`` squared distances between two batched point sets.

        The Gram form rather than a broadcast difference: it is one BLAS call instead of
        an (B, na, nb, 3) temporary, and for the magnitudes involved here (coordinates of
        order 20 A, distances of order 3 A) the cancellation costs ~1e-14 A, four orders
        below the tolerance `selftest` enforces.
        """
        g = A @ np.swapaxes(B, 1, 2)
        # Accumulated in place: the naive expression allocates three (B, na, nb)
        # temporaries, and at these sizes the model is memory-bound, not FLOP-bound.
        np.multiply(g, -2.0, out=g)
        g += (A * A).sum(2)[:, :, None]
        g += (B * B).sum(2)[:, None, :]
        return g

    def _steric(self, atoms: np.ndarray) -> np.ndarray:
        # Every operation is in place on the one (B, A, A) buffer. Gathering the ~830
        # listed pairs out of it first was measured 40% SLOWER than just running the
        # arithmetic over the full matrix: fancy indexing costs more per element than the
        # dense passes it saves. Non-listed entries have limit 0, and no distance is
        # negative, so their overlap is identically 0.
        d = self._sqdist(atoms, atoms)
        np.maximum(d, 0.0, out=d)
        np.sqrt(d, out=d)
        np.subtract(self._st_lim_mat, d, out=d)
        np.maximum(d, 0.0, out=d)
        return np.einsum("bij,bij->b", d, d)

    @staticmethod
    def _switch(d: np.ndarray, d0: float, dc: float) -> np.ndarray:
        s = np.zeros_like(d)
        s[d <= d0] = 1.0
        mid = (d > d0) & (d < dc)
        s[mid] = 0.5 * (1.0 + np.cos(math.pi * (d[mid] - d0) / (dc - d0)))
        return s

    def _rama(self, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
        """(B,) Ramachandran sum, matching `energy_terms.rama_penalty` term by term."""
        pd = np.degrees(phi)
        sd = np.degrees(psi)
        dphi = ((pd[:, :, None] - self._rama_pc[None, None, :] + 180.0) % 360.0) - 180.0
        dpsi = ((sd[:, :, None] - self._rama_sc[None, None, :] + 180.0) % 360.0) - 180.0
        d2 = dphi * dphi + dpsi * dpsi
        score = (self._rama_W[None] * np.exp(-d2 / self._rama_two_sig2[None, None, :])
                 ).sum(2)
        e = 1.0 - score
        if self._is_pro.any():
            p = pd[:, self._is_pro]
            e[:, self._is_pro] += (0.03 * np.maximum(0.0, -90.0 - p)
                                   + 0.03 * np.maximum(0.0, p + 50.0))
        return e.sum(1) + self._gly_shift

    def _hbonds(self, coords: Dict[str, np.ndarray]
                ) -> Tuple[np.ndarray, np.ndarray, List[Tuple[Tuple[int, int], ...]]]:
        """Batched DSSP matrix, then the same stable greedy match as the scalar path."""
        N, C, O = coords["N"], coords["C"], coords["O"]
        B, n = N.shape[0], N.shape[1]
        H = np.full((B, n, 3), np.nan)
        if n >= 2:
            d = C[:, :-1] - O[:, :-1]
            nd = np.linalg.norm(d, axis=2)
            ok_h = nd > 1e-6
            Hn = N[:, 1:] + np.where(ok_h[..., None], d / np.where(ok_h, nd, 1.0)[..., None],
                                     np.nan)
            H[:, 1:] = Hn
        valid = np.isfinite(H).all(axis=2)

        with np.errstate(divide="ignore", invalid="ignore"):
            # Index convention follows `protein_geometry.dssp_energy_matrix` exactly:
            # row = donor residue (its N and H), column = acceptor residue (its C and O).
            dON = np.sqrt(np.maximum(self._sqdist(N, O), 0.0))
            dCH = np.sqrt(np.maximum(self._sqdist(H, C), 0.0))
            dOH = np.sqrt(np.maximum(self._sqdist(H, O), 0.0))
            dCN = np.sqrt(np.maximum(self._sqdist(N, C), 0.0))
            E = _HB_PREFACTOR * (1.0 / dON + 1.0 / dCH - 1.0 / dOH - 1.0 / dCN)
        E = np.maximum(E, geo.HB_E_FLOOR)

        ok = valid[:, :, None] & self._hb_sep_ok[None]
        with np.errstate(invalid="ignore"):
            ok &= (dON > geo.HB_MIN_ON) & (dCH > 0.5) & (dOH > 0.5) & (dCN > 0.5)
            ok &= np.isfinite(E) & (E < -0.5)

        local = np.zeros(B)
        lr = np.zeros(B)
        all_pairs: List[Tuple[Tuple[int, int], ...]] = [()] * B
        rows, dons, accs = np.nonzero(ok)          # row-major, same order as np.where
        if rows.size == 0:
            return local, lr, all_pairs
        en_all = E[rows, dons, accs]
        # One global sort that is stable *within* each structure by energy reproduces the
        # scalar path's per-structure `argsort(..., kind="stable")` exactly: `np.nonzero`
        # already emits row-major order, and the arange key preserves it among ties (which
        # HB_E_FLOOR manufactures at -4.0, so the tie-break is not hypothetical).
        order = np.lexsort((np.arange(rows.size), en_all, rows))
        rows_s = rows[order]
        dd = dons[order].tolist()
        aa = accs[order].tolist()
        en = en_all[order].tolist()
        bounds = np.searchsorted(rows_s, np.arange(B + 1))
        sep_cut = et.HB_LONGRANGE_SEP
        for b in range(B):
            lo, hi = int(bounds[b]), int(bounds[b + 1])
            if lo == hi:
                continue
            # Bitmasks rather than boolean arrays: the loop body runs once per admissible
            # pair over the whole batch, and two small np.zeros calls per structure cost
            # more than the matching itself.
            du = au = 0
            matched: List[Tuple[int, int]] = []
            lo_s = lr_s = 0.0
            for k in range(lo, hi):
                i, j = dd[k], aa[k]
                bi, bj = 1 << i, 1 << j
                if (du & bi) or (au & bj):
                    continue
                du |= bi
                au |= bj
                e = en[k] + 1.0                   # desolvation_cost
                if e >= 0.0:
                    continue
                matched.append((i, j))
                if abs(i - j) < sep_cut:
                    lo_s += e
                else:
                    lr_s += e
            local[b] = lo_s
            lr[b] = lr_s
            all_pairs[b] = tuple(matched)
        return local, lr, all_pairs

    # -- public ------------------------------------------------------------
    def components(self, phi: np.ndarray, psi: np.ndarray,
                   coords: Optional[Dict[str, np.ndarray]] = None
                   ) -> Dict[str, np.ndarray]:
        """Unweighted term values, each a ``(B,)`` array, for ``(B, n)`` torsions."""
        phi = np.atleast_2d(np.asarray(phi, dtype=float))
        psi = np.atleast_2d(np.asarray(psi, dtype=float))
        if coords is None:
            coords = geo.build_backbone_batch(phi, psi)
        CA, CB = coords["CA"], coords["CB"]
        B = CA.shape[0]
        di, dj = self._di, self._dj
        d_cb = np.linalg.norm(CB[:, di] - CB[:, dj], axis=2)

        atoms = np.concatenate([coords[k] for k in _ATOM_ORDER], axis=1)

        contact = self._switch(d_cb[:, self._m3], 4.5, 8.5) @ self._mj3
        solvation = -(self._switch(d_cb, 6.0, 10.0) @ self._sol_w)

        if self._el_qq.size:
            d = d_cb[:, self._el_mask]
            el = ((et.COULOMB / et.DIELECTRIC)
                  * (self._el_qq[None] / np.maximum(d, 2.0))
                  * np.exp(-d / 8.0)).sum(1)
        else:
            el = np.zeros(B)

        if self._arom_i.size:
            da = np.linalg.norm(CB[:, self._arom_i] - CB[:, self._arom_j], axis=2)
            arom = -np.exp(-((da - et.AROM_CB_DIST) / et.AROM_CB_WIDTH) ** 2).sum(1)
        else:
            arom = np.zeros(B)

        hb_local, hb_lr, pairs = self._hbonds(coords)
        coop_h = np.zeros(B)
        coop_s = np.zeros(B)
        for b, p in enumerate(pairs):
            if p:                       # both terms are 0 on an empty pair list
                coop_h[b] = et.coop_helix_term(p)
                coop_s[b] = et.coop_sheet_term(p)

        c = CA - CA.mean(1, keepdims=True)
        rg = np.sqrt((c * c).sum(2).mean(1))
        compact = np.maximum(0.0, rg - self._rg_target) ** 2

        return {
            "steric": self._steric(atoms),
            "contact": contact,
            "hbond_local": hb_local,
            "hbond_longrange": hb_lr,
            "coop_helix": coop_h,
            "coop_sheet": coop_s,
            "solvation": solvation,
            "electrostatic": el,
            "aromatic": arom,
            "torsion": self._rama(phi, psi),
            "compactness": compact,
        }

    def total(self, comp: Dict[str, np.ndarray],
              weights: Optional[Dict[str, float]] = None) -> np.ndarray:
        w = self.weights if weights is None else weights
        out = np.zeros_like(comp["steric"])
        for k in et.TERM_NAMES:
            out = out + w.get(k, 0.0) * comp[k]
        return out

    def energy(self, phi: np.ndarray, psi: np.ndarray,
               coords: Optional[Dict[str, np.ndarray]] = None) -> np.ndarray:
        """``(B,)`` weighted Legacy energy for ``(B, n)`` phi/psi in radians."""
        return self.total(self.components(phi, psi, coords))

    def component_matrix(self, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
        """``(B, 11)`` unweighted components in `energy_terms.TERM_NAMES` order."""
        comp = self.components(phi, psi)
        return np.column_stack([comp[k] for k in et.TERM_NAMES])


def batched_energy(sequence: str, phi: np.ndarray, psi: np.ndarray,
                   weights: Optional[Dict[str, float]] = None) -> np.ndarray:
    """One-shot convenience wrapper; build a `BatchedLegacy` if calling repeatedly."""
    return BatchedLegacy(sequence, weights).energy(phi, psi)


# ==========================================================================
# Weight vectors (for fitting)
# ==========================================================================
def vector_from_weights(weights: Dict[str, float],
                        names: Sequence[str] = et.TERM_NAMES) -> np.ndarray:
    return np.array([float(weights.get(k, 0.0)) for k in names])


def weights_from_vector(vec: Iterable[float],
                        names: Sequence[str] = et.TERM_NAMES,
                        base: Optional[Dict[str, float]] = None) -> Dict[str, float]:
    w = dict(et.DEFAULT_WEIGHTS if base is None else base)
    for k, v in zip(names, vec):
        w[k] = float(v)
    return w


# ==========================================================================
# Refinement
# ==========================================================================
def _clip_box(phi, psi, phi0, psi0, box):
    return (np.clip(phi, phi0 - box, phi0 + box),
            np.clip(psi, psi0 - box, psi0 + box))


def refine_angles(sequence: str, phi0: np.ndarray, psi0: np.ndarray,
                  weights: Optional[Dict[str, float]] = None,
                  box_deg: float = DEFAULT_BOX_DEG,
                  steps: int = DEFAULT_STEPS,
                  pop: int = 64, seed: int = 0,
                  polish: bool = True,
                  model: Optional[BatchedLegacy] = None) -> Dict[str, object]:
    """Bounded local minimisation of the Legacy energy over continuous (phi, psi).

    Deterministic: the proposal stream is a fixed `default_rng(seed)` sequence and the
    accept rule is a strict improvement test, so the same inputs give the same output.

    The search is a batched pattern search rather than a gradient descent because the
    objective has integer-valued pieces (`coop_helix`, `coop_sheet`) and a hard admission
    cutoff in `hbond_*`; a finite-difference gradient is exactly zero across those
    plateaus. `polish` adds an L-BFGS-B pass afterwards, which does help on the smooth
    remainder once the discrete pieces have settled.

    Returns ``{phi, psi, coords, CA, energy, energy0, evals}`` -- angles in radians.
    """
    m = model or BatchedLegacy(sequence, weights)
    phi0 = np.asarray(phi0, dtype=float).ravel()
    psi0 = np.asarray(psi0, dtype=float).ravel()
    box = math.radians(float(box_deg))
    rng = np.random.default_rng(seed)
    n = m.n

    cur_phi, cur_psi = phi0.copy(), psi0.copy()
    e0 = float(m.energy(cur_phi[None], cur_psi[None])[0])
    best = e0
    evals = 1

    sigma = box * 0.6
    sigma_min = box * 0.02
    for _ in range(steps):
        d = rng.normal(0.0, sigma, size=(pop, 2 * n))
        cand_phi, cand_psi = _clip_box(cur_phi[None] + d[:, :n],
                                       cur_psi[None] + d[:, n:], phi0, psi0, box)
        e = m.energy(cand_phi, cand_psi)
        evals += pop
        k = int(np.argmin(e))
        if e[k] < best - 1e-12:
            best = float(e[k])
            cur_phi, cur_psi = cand_phi[k].copy(), cand_psi[k].copy()
        else:
            sigma *= 0.75
            if sigma < sigma_min:
                break

    if polish:
        try:
            from scipy.optimize import minimize
        except ImportError:
            minimize = None
        if minimize is not None:
            x0 = np.concatenate([cur_phi, cur_psi])
            lo = np.concatenate([phi0 - box, psi0 - box])
            hi = np.concatenate([phi0 + box, psi0 + box])
            counter = {"n": 0}
            h = 1e-4
            eye = np.eye(2 * n) * h

            def fg(x):
                # The point and its 4n central-difference neighbours in ONE batched call.
                # Letting L-BFGS-B take the differences itself would issue 4n separate
                # B=1 evaluations, where the fixed numpy overhead per call dominates.
                X = np.vstack([x[None], x[None] + eye, x[None] - eye])
                e = m.energy(X[:, :n], X[:, n:])
                counter["n"] += len(X)
                g = (e[1:2 * n + 1] - e[2 * n + 1:]) / (2.0 * h)
                return float(e[0]), g

            res = minimize(fg, x0, method="L-BFGS-B", jac=True,
                           bounds=list(zip(lo, hi)),
                           options={"maxiter": 100, "ftol": 1e-12, "gtol": 1e-8})
            evals += counter["n"]
            if float(res.fun) < best:
                best = float(res.fun)
                cur_phi, cur_psi = np.clip(res.x[:n], lo[:n], hi[:n]), \
                    np.clip(res.x[n:], lo[n:], hi[n:])

    coords = geo.build_backbone_batch(cur_phi[None], cur_psi[None])
    return {"phi": cur_phi, "psi": cur_psi,
            "coords": {k: v[0] for k, v in coords.items()},
            "CA": coords["CA"][0], "energy": best, "energy0": e0, "evals": evals}


def refine(sequence: str, rep, states, weights: Optional[Dict[str, float]] = None,
           box_deg: float = DEFAULT_BOX_DEG, steps: int = DEFAULT_STEPS,
           pop: int = 64, seed: int = 0, polish: bool = True,
           model: Optional[BatchedLegacy] = None) -> Dict[str, object]:
    """Refine the structure encoded by discrete `states` under `rep`.

    `states` is either an ``(n,)`` array of per-residue state indices or a bitstring the
    representation can decode.
    """
    if isinstance(states, str):
        phi0, psi0, _, _ = rep.build_all(states)
        phi0 = np.asarray(phi0, dtype=float)
        psi0 = np.asarray(psi0, dtype=float)
    else:
        s = np.asarray(states, dtype=int).ravel()
        rows = np.arange(len(s))
        phi0 = rep._phi[rows, s]
        psi0 = rep._psi[rows, s]
    return refine_angles(sequence, phi0, psi0, weights=weights, box_deg=box_deg,
                         steps=steps, pop=pop, seed=seed, polish=polish, model=model)


# ==========================================================================
# Self-test
# ==========================================================================
def selftest(sequence: str = "RLKWVRIWRRGDYE", B: int = 24, seed: int = 0,
             tol: float = 1e-9, verbose: bool = True) -> float:
    """Assert the batched terms match `energy_terms.energy_components` to `tol`.

    Random torsions rather than library states, so the comparison exercises continuous
    angles -- the regime `refine` runs in and the one where `rama_penalty`'s memoisation
    stops hiding a mismatch.
    """
    rng = np.random.default_rng(seed)
    n = len(sequence)
    phi = rng.uniform(-math.pi, math.pi, (B, n))
    psi = rng.uniform(-math.pi, math.pi, (B, n))
    coords = geo.build_backbone_batch(phi, psi)
    m = BatchedLegacy(sequence)
    fast = m.components(phi, psi, coords)
    worst = 0.0
    for b in range(B):
        co = {k: coords[k][b] for k in coords}
        ref = et.energy_components(sequence, co, phi[b], psi[b])
        for k in et.TERM_NAMES:
            d = abs(float(fast[k][b]) - ref[k])
            scale = max(1.0, abs(ref[k]))
            worst = max(worst, d / scale)
            if d / scale > tol:
                raise AssertionError(
                    f"term {k} mismatch at b={b}: batched {fast[k][b]!r} "
                    f"vs scalar {ref[k]!r}")
    if verbose:
        print(f"selftest ok: worst relative term error {worst:.3e} "
              f"over {B} structures x {len(et.TERM_NAMES)} terms")
    return worst


if __name__ == "__main__":
    selftest()
    selftest("GYDPETGTWG", B=32, seed=3)
    selftest("PGWALCDEFK", B=32, seed=5)
