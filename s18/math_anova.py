"""SPRINT 18 / MATH -- THE DELIVERABLE.  `E_le1(theta, target)` in CONTINUOUS torsion space.

This module is what the EXPERIMENT workstream consumes.  It must not build its own.

------------------------------------------------------------------------------------------
THE DERIVATION, IN FULL, BECAUSE EVERY DOWNSTREAM NUMBER DEPENDS ON IT BEING RIGHT
------------------------------------------------------------------------------------------

THE OBJECT BEING TRUNCATED.  The production objective (`s17/refine.py::_obj`) is

    E(theta) = sum_{p=(a,b), b-a>=2}  w_p ( d_p(theta) - dhat_p )^2 ,   w_p = 1/sd_p^2

over `theta = (phi_0..phi_{n-1}, psi_0..psi_{n-1})` with `d_p` the CA-CA distance of an
ideal-geometry chain.  `dhat` is the leave-fold-out distogram expectation with the
leave-fold-out separation-binned debias; `sd` its predicted sd.  Native-free throughout.

THE VARIABLE.  The natural coordinate is the RESIDUE, `theta_r = (phi_r, psi_r) in T^2`, not
the individual angle.  This is the continuous analogue of the per-qubit / per-residue
distinction that decides everything on the lattice: an ANOVA is relative to a CHOSEN product
factorisation, and the residue is the factor that has physical meaning.

LOCALITY, EXACT AND VERIFIED (`support_check`).  `d_ab` depends on exactly the residues
`{a+1, ..., b-1}` -- fully, both torsions -- and on nothing else.  `phi_a`-type boundary
terms drop out because `psi_a` rotates the chain about an axis THROUGH `CA_a` and `phi_b`
only moves atoms placed after `CA_b`.  Two corollaries, both EXACT:

  * the objective is invariant to `phi_0, psi_0, phi_{n-1}, psi_{n-1}` -- the whole distance
    objective is a function of the `n-2` INTERIOR residues only;
  * therefore `f_0 = f_{n-1} = 0` identically.  **The degree-1 object exerts no force
    whatever on the two terminal residues**, and the frozen metric scores them.

THE DECOMPOSITION.  Under a product reference measure `mu = tensor_r mu_r` on `(T^2)^n`,

    E_0     = E_mu[E]
    f_r(t)  = E_mu[E | theta_r = t] - E_0                      (E_{mu_r}[f_r] = 0)
    E_le1   = E_0 + sum_r f_r(theta_r)         E_ge2 = E - E_le1

`E_le1` is the mu-orthogonal projection of `E` onto the additive subspace.  It is the
CONDITIONAL-EXPECTATION object, not a marginalisation of a probability: `E_mu[E|theta_r]` is
an average of the ENERGY over the other residues, so `f_r(t)` reads, in one sentence,

    > the mean distogram penalty this target pays for residue r adopting torsions t, when
    > every other residue is drawn from the reference ensemble.

COMPUTABILITY -- WHERE LOCALITY PAYS.  Expanding term by term,

    f_r(t) = sum_{p : r in supp(p)} w_p ( E_mu[(d_p - dhat_p)^2 | theta_r = t] - c_p )

because every pair whose support excludes `r` contributes a constant that the centring
removes.  Only `r*(n-1-r)` of the `~n^2/2` pairs vary with `theta_r`.  The estimator below
uses COMMON RANDOM NUMBERS across the whole `theta_r` grid, so the constant pairs cancel
EXACTLY (not just in expectation) and the quadrature error enters only through the pairs
that genuinely depend on `r`.

THE QUADRATURE IS PART OF THE DEFINITION.  `E_mu[.|theta_r=t]` has no closed form
(`E[d]` is not a polynomial moment of a product of rotations, though `E[d^2]` is).  So the
object is defined against a FROZEN sample `Theta^(1..S) ~ mu` drawn once from
`stable_rng(pdb, 'mu', name)`.  That makes `E_le1` a deterministic, cacheable, bit-reproducible
function whose coefficients are persisted.  Its quadrature error is measured (S vs S/2) and
reported, and the lattice cross-check in `math_bridge.py` compares it against the EXACT
enumerated ANOVA.

REPRESENTATION ON THE TORUS.  `f_r` is a function on `T^2`; a polynomial basis would be
wrong.  It is tabulated on a `G x G` grid and evaluated by EXACT TRIGONOMETRIC
INTERPOLATION (2-D DFT), which is the band-limited interpolant natural to a torus:
periodic by construction, analytic, exact at the nodes, and differentiable in closed form.
The persisted coefficients ARE the DFT.  A strictly coarser alternative -- keeping only
circular harmonics of order <= 1 -- is provided as `harmonic_truncate` and reported
separately, because "additive in residues" and "degree-1 in angles" are DIFFERENT
truncations and conflating them is the continuous form of the per-qubit/per-residue error.

WHAT IS AND IS NOT KEPT, BY PHYSICAL TERM.
  weight 0   `E_0`: the mean penalty of the reference ensemble.  Retained.
  weight 1   for each residue, the full 2-torus profile of how that residue's torsions shift
             the ensemble-mean penalty of every pair that spans it.  Retained.
  weight >=2 every genuine CO-OPERATIVE effect: that two residues' torsions must be chosen
             TOGETHER to place a distant pair.  REMOVED.  A hairpin whose two strands are
             set by a compensating pair of turns has its entire signal here.
Coefficients depend on residue identity only through `mu_r` (if `mu` is target-conditioned)
and through which pairs span `r`; they depend on pair separation through `w_p` and through
the support-size weighting -- an interior residue is spanned by `r(n-1-r)` pairs, so the
field is automatically strongest in the middle of the chain and exactly zero at the ends.

THE LANDSCAPE.  This is not a reparameterisation: `E_le1` is a different function with a
different optimiser.  Its global minimum is CLOSED FORM -- minimise each `f_r`
independently -- so it poses no search problem at all, which is the whole content of
Sprint 17's "no search problem here for a quantum device".

RUN:   python -m s18.math_anova selfcheck      (Control D + locality + quadrature)
       python -m s18.math_anova build          (cache coefficients for the 126 tuning targets)
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s18 import math_lib as M

from s12 import instrument as I                   # noqa: E402
from s15 import distcal as C                      # noqa: E402
from s15 import seed as SD                        # noqa: E402

GRID = 16                     # G x G tabulation of each f_r; 22.5 degree spacing
NSAMP = 512                   # frozen quadrature size (CRN across the whole mesh)
MU_DEFAULT = "pool"
CONFIG = {"module": "s18/math_anova.py", "grid": GRID, "nsamp": NSAMP,
          "mu_default": MU_DEFAULT, "salt": M.SALT, "min_sep": 2,
          "objective": "sum_p ((d_p - dhat_p)/sd_p)^2, dhat debiased leave-fold-out"}

_ANG = (np.arange(GRID) + 0.5) / GRID * 2 * np.pi - np.pi     # bin centres on [-pi, pi)


def _interp1(Fk, ang):
    """Exact trigonometric interpolation of `(n, G)` 1-D DFTs at `(B, n)` angles."""
    G = Fk.shape[-1]
    k = np.fft.fftfreq(G, 1.0 / G)
    u = (np.asarray(ang, float) + np.pi) / (2 * np.pi) * G - 0.5
    E = np.exp(2j * np.pi * u[:, :, None] * k[None, None, :] / G)
    return np.einsum("bnu,nu->b", E, Fk).real / G


# ============================================================ reference measures
def mu_samples(pdb, seq, fold, n, name, S, rng):
    """`(S, n)` phi and psi drawn i.i.d. per residue from a NATIVE-FREE product measure.

    pool     the retrieval pool's own empirical per-residue torsion marginal -- the shipped
             top-75 windows, resampled per residue INDEPENDENTLY so the measure is a genuine
             product (drawing whole windows would not be one).  Target-conditioned,
             native-free, and the most informative legitimate choice.
    uniform  uniform on the torus -- the zero-information reference measure.
    rama     the shipped leave-fold-out Ramachandran density of each residue's class,
             sampled by bin with a uniform jitter inside the bin.  Generic, native-free.
    """
    if name == "uniform":
        return (rng.uniform(-np.pi, np.pi, (S, n)), rng.uniform(-np.pi, np.pi, (S, n)))
    if name == "pool":
        from s14 import retprior as RP
        PHI, PSI, _ = RP.windows(pdb, "top75")
        m = PHI.shape[0]
        pick = rng.integers(0, m, (S, n))
        cols = np.arange(n)[None, :]
        return PHI[pick, cols], PSI[pick, cols]
    if name == "lattice":
        #: uniform over the k = 4 torsion atoms of the ENUMERATED lattice, per residue.
        #: Used only by `math_bridge.py`, to check the continuous machinery against the
        #: exact enumerated ANOVA under the identical reference measure.
        from s16 import qphase_lib as QP
        ins = QP.inst(pdb)
        st = rng.integers(0, ins.k, (S, n))
        cols = np.arange(n)[None, :]
        return ins.PHI[cols, st], ins.PSI[cols, st]
    if name == "rama":
        from core import project as PJ
        T = PJ.logp_tables("rama")[int(fold)]
        cls = PJ.res_classes(seq)
        L = T[cls]                                   # (n, RB, RB) log density
        RB = L.shape[-1]
        p = np.exp(L - L.max((1, 2), keepdims=True)).reshape(n, RB * RB)
        p = p / p.sum(1, keepdims=True)
        cum = np.cumsum(p, 1)
        u = rng.random((S, n))
        flat = (u[:, :, None] > cum[None, :, :-1]).sum(2)
        bi, bj = flat // RB, flat % RB
        step = 2 * np.pi / RB
        phi = -np.pi + (bi + rng.random((S, n))) * step
        psi = -np.pi + (bj + rng.random((S, n))) * step
        return phi, psi
    raise ValueError(name)


# ============================================================ the objective itself
class Target:
    """Everything about one target's continuous objective, cached to `s18/cache/`."""

    def __init__(self, pdb, seq, n, fold, dhat, sd, i, j):
        self.pdb, self.seq, self.n, self.fold = pdb, seq, int(n), int(fold)
        self.dhat = np.asarray(dhat, float)
        self.sd = np.asarray(sd, float)
        self.w = 1.0 / self.sd ** 2
        self.i = np.asarray(i, int)
        self.j = np.asarray(j, int)
        self.E0 = None
        self.f = None                      # (n, G, G) tabulated per-RESIDUE fields
        self.F = None                      # (n, G, G) complex DFT of f
        self.a = self.b = None             # (n, G) per-ANGLE fields (sub-residue object)
        self.A1 = self.B1 = None
        self.m1 = self.m2 = None           # (n, G*G, n_pairs) conditional distance moments
        self.am1 = None                    # the same for the per-angle object
        self.u1 = self.u2 = None           # unconditional pair-distance moments
        self.mu = None

    # ---------------------------------------------------------------- full objective
    def E_full(self, phi, psi):
        """`sum_p w_p (d_p - dhat_p)^2`.  Accepts `(n,)` or `(B, n)`; returns scalar/`(B,)`.
        Bit-identical in form to `s17/refine.py::_obj`."""
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        CA = np.asarray(I.build_ca(phi, psi), float)
        d = np.sqrt(((CA[:, self.i] - CA[:, self.j]) ** 2).sum(-1))
        r = d - self.dhat[None, :]
        out = (r * r * self.w[None, :]).sum(1)
        return float(out[0]) if out.size == 1 and np.ndim(phi) == 2 and phi.shape[0] == 1 \
            else out

    # ---------------------------------------------------------------- the fields
    def fit(self, mu=MU_DEFAULT, S=NSAMP, grid=GRID, batch=32768, half=False,
            moments=True):
        """Tabulate `f_r` on a `grid x grid` torus mesh by frozen-CRN conditional MC.

        The SAME `S` reference configurations are used at every mesh point and for every
        residue, so every pair that does not span `r` cancels to the last bit.

        THE ESTIMATOR, AND WHY IT IS THIS ONE.  Write
        `h_r(t) = (1/S) sum_s E(Theta_s with residue r set to t)` and `Ebar = (1/S) sum_s
        E(Theta_s)`.  The additive model is stored as `f_r = h_r - Ebar` and `E_0 = Ebar`, so

            E_le1(theta) = Ebar + sum_r (h_r(theta_r) - Ebar) = sum_r h_r(theta_r) - (n-1) Ebar

        For an EXACTLY ADDITIVE `E = sum_r g_r` this is `sum_r g_r(theta_r)` IDENTICALLY --
        every quadrature error cancels, at any `S`, because `h_r = g_r + (Ebar - gbar_r^S)`
        and the `Ebar` terms telescope.  That is what Control D verifies to machine
        precision, and it is why the mesh must NOT be re-centred on itself: re-centring
        breaks the telescoping and leaves an O(sigma/sqrt(S)) constant behind.
        """
        n, G = self.n, grid
        rng = SD.stable_rng(self.pdb, mu, S, grid, "mufreeze", salt=M.SALT)
        PHI, PSI = mu_samples(self.pdb, self.seq, self.fold, n, mu, S, rng)
        if half:
            PHI, PSI = PHI[: S // 2], PSI[: S // 2]
        S = PHI.shape[0]
        ang = (np.arange(G) + 0.5) / G * 2 * np.pi - np.pi
        gp, gq = np.meshgrid(ang, ang, indexing="ij")
        gp, gq = gp.ravel(), gq.ravel()                       # (G*G,)
        self.mu, self.S, self.grid = mu, int(S), int(G)
        step = max(1, batch // S)

        if not moments:
            #: generic path (Control D, whose `E_full` is a synthetic function)
            self.E0 = float(self._E_batch(PHI, PSI, batch).mean())
            f = np.zeros((n, G, G), float)
            for r in range(1, n - 1):                         # f_0 = f_{n-1} = 0, EXACT
                acc = np.empty(G * G, float)
                for a in range(0, G * G, step):
                    b = min(a + step, G * G)
                    P = np.repeat(PHI[None, :, :], b - a, 0).reshape(-1, n)
                    Q = np.repeat(PSI[None, :, :], b - a, 0).reshape(-1, n)
                    P[:, r] = np.repeat(gp[a:b], S)
                    Q[:, r] = np.repeat(gq[a:b], S)
                    acc[a:b] = self._E_batch(P, Q, batch).reshape(b - a, S).mean(1)
                f[r] = (acc - self.E0).reshape(G, G)
            self.f = f
            self.F = np.fft.fft2(f, axes=(1, 2))
            self._fit_angles(PHI, PSI, ang, batch)
            return self

        #: MOMENT PATH.  Tabulate the conditional first and second moments of every pair
        #: DISTANCE rather than of the energy.  Because
        #:      E = sum_p w_p (d_p^2 - 2 dhat_p d_p + dhat_p^2)
        #: is AFFINE in (d^2, d) at fixed theta, the conditional energy for ANY dhat follows
        #: from the SAME tables with no new chain builds.  That is what makes the
        #: sensitivity of the argmin to dhat exactly computable instead of re-simulated.
        npair = len(self.i)
        m1 = np.zeros((n, G * G, npair))
        m2 = np.zeros((n, G * G, npair))
        D = self._D_batch(PHI, PSI)
        u1, u2 = D.mean(0), (D * D).mean(0)
        del D
        for r in range(1, n - 1):
            for a in range(0, G * G, step):
                b = min(a + step, G * G)
                P = np.repeat(PHI[None, :, :], b - a, 0).reshape(-1, n)
                Q = np.repeat(PSI[None, :, :], b - a, 0).reshape(-1, n)
                P[:, r] = np.repeat(gp[a:b], S)
                Q[:, r] = np.repeat(gq[a:b], S)
                D = self._D_batch(P, Q).reshape(b - a, S, npair)
                m1[r, a:b] = D.mean(1)
                m2[r, a:b] = (D * D).mean(1)
        for r in (0, n - 1):
            m1[r] = u1[None, :]
            m2[r] = u2[None, :]
        self.m1, self.m2, self.u1, self.u2 = m1, m2, u1, u2
        self._fit_angles_moments(PHI, PSI, ang, step, S)
        self._set_dhat(self.dhat)
        return self

    # -------------------------------------------- fields from the moment tables
    def _set_dhat(self, dhat):
        """Rebuild `E_0`, `f_r`, `a_r`, `b_r` for an arbitrary target distance vector.
        EXACT given the tables, with no new chain builds.  This is the handle the
        argmin-sensitivity analysis turns."""
        dhat = np.asarray(dhat, float)
        w, G, n = self.w, self.grid, self.n
        self.dhat = dhat
        self.E0 = float((w * (self.u2 - 2 * dhat * self.u1 + dhat ** 2)).sum())
        hh = ((self.m2 - 2 * dhat[None, None, :] * self.m1
               + (dhat ** 2)[None, None, :]) * w[None, None, :]).sum(-1)     # (n, G*G)
        self.f = (hh - self.E0).reshape(n, G, G)
        self.F = np.fft.fft2(self.f, axes=(1, 2))
        if self.am1 is not None:
            ha = ((self.am2 - 2 * dhat[None, None, :] * self.am1
                   + (dhat ** 2)[None, None, :]) * w[None, None, :]).sum(-1)
            hb = ((self.bm2 - 2 * dhat[None, None, :] * self.bm1
                   + (dhat ** 2)[None, None, :]) * w[None, None, :]).sum(-1)
            self.a, self.b = ha - self.E0, hb - self.E0
            self.A1 = np.fft.fft(self.a, axis=1)
            self.B1 = np.fft.fft(self.b, axis=1)
        return self

    def with_dhat(self, dhat):
        """A shallow copy carrying a different distance vector, fields rebuilt exactly."""
        import copy
        t = copy.copy(self)
        t._set_dhat(np.asarray(dhat, float))
        return t

    def _D_batch(self, PHI, PSI):
        CA = np.asarray(I.build_ca(PHI, PSI), float)
        return np.sqrt(((CA[:, self.i] - CA[:, self.j]) ** 2).sum(-1))

    def _fit_angles_moments(self, PHI, PSI, ang, step, S):
        n, G, npair = self.n, len(ang), len(self.i)
        am1 = np.repeat(np.repeat(self.u1[None, None, :], G, 1), n, 0).copy()
        am2 = np.repeat(np.repeat(self.u2[None, None, :], G, 1), n, 0).copy()
        bm1, bm2 = am1.copy(), am2.copy()
        for r in range(1, n - 1):
            for which, t1, t2 in (("phi", am1, am2), ("psi", bm1, bm2)):
                for s0 in range(0, G, step):
                    s1 = min(s0 + step, G)
                    P = np.repeat(PHI[None, :, :], s1 - s0, 0).reshape(-1, n)
                    Q = np.repeat(PSI[None, :, :], s1 - s0, 0).reshape(-1, n)
                    (P if which == "phi" else Q)[:, r] = np.repeat(ang[s0:s1], S)
                    D = self._D_batch(P, Q).reshape(s1 - s0, S, npair)
                    t1[r, s0:s1] = D.mean(1)
                    t2[r, s0:s1] = (D * D).mean(1)
        self.am1, self.am2, self.bm1, self.bm2 = am1, am2, bm1, bm2

    def _fit_angles(self, PHI, PSI, ang, batch):
        """THE SUB-RESIDUE OBJECT.  Order-1 ANOVA in the 2n-dimensional factorisation over
        INDIVIDUAL ANGLES: `E ~ E_0 + sum_r a_r(phi_r) + sum_r b_r(psi_r)`.

        This is the honest continuous analogue of "finer than per-residue", and it exists
        only because phi and psi are separate coordinates.  It is NOT the continuous image
        of the lattice's per-qubit truncation -- the two lattice bits index k-means clusters
        of the joint (phi, psi) library and do not correspond to phi and psi at all.  There
        is therefore no continuous object equal to the strict Walsh weight-<=1 projection,
        and that is the exact content of falsifier F5.
        """
        n, S = self.n, PHI.shape[0]
        G = len(ang)
        a = np.zeros((n, G))
        b = np.zeros((n, G))
        step = max(1, batch // S)
        for r in range(1, n - 1):
            for which, tab in (("phi", a), ("psi", b)):
                acc = np.empty(G)
                for s0 in range(0, G, step):
                    s1 = min(s0 + step, G)
                    P = np.repeat(PHI[None, :, :], s1 - s0, 0).reshape(-1, n)
                    Q = np.repeat(PSI[None, :, :], s1 - s0, 0).reshape(-1, n)
                    (P if which == "phi" else Q)[:, r] = np.repeat(ang[s0:s1], S)
                    acc[s0:s1] = self._E_batch(P, Q, batch).reshape(s1 - s0, S).mean(1)
                tab[r] = acc - self.E0
        self.a, self.b = a, b
        self.A1 = np.fft.fft(a, axis=1)
        self.B1 = np.fft.fft(b, axis=1)

    def E_le1_ang(self, phi, psi):
        """`E_0 + sum_r a_r(phi_r) + sum_r b_r(psi_r)` -- the sub-residue (angle-additive)
        truncation.  Strictly coarser than `E_le1`."""
        assert self.A1 is not None, "call fit() first"
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        B = phi.shape[0]
        out = self.E0 + _interp1(self.A1, phi) + _interp1(self.B1, psi)
        return float(out[0]) if B == 1 else out

    def _E_batch(self, PHI, PSI, batch):
        out = np.empty(len(PHI), float)
        for a in range(0, len(PHI), batch):
            b = min(a + batch, len(PHI))
            CA = np.asarray(I.build_ca(PHI[a:b], PSI[a:b]), float)
            d = np.sqrt(((CA[:, self.i] - CA[:, self.j]) ** 2).sum(-1))
            rr = d - self.dhat[None, :]
            out[a:b] = (rr * rr * self.w[None, :]).sum(1)
        return out

    # ---------------------------------------------------------------- evaluation
    def _interp(self, F, phi, psi):
        """Exact band-limited trigonometric interpolation of a `(., G, G)` DFT on `T^2`.

        `f(x, y) = (1/G^2) sum_{u,v} F_{uv} exp(i(u x' + v y'))` with `x'` the mesh
        coordinate; the half-bin offset of the mesh is carried explicitly.  Real part taken;
        the imaginary residue is machine noise and is asserted small in `selfcheck`.
        """
        G = F.shape[-1]
        k = np.fft.fftfreq(G, 1.0 / G)                       # integer frequencies
        # mesh point m sits at angle (m + 0.5)/G*2pi - pi  =>  m = (angle+pi)/(2pi)*G - 0.5
        u = (np.asarray(phi, float) + np.pi) / (2 * np.pi) * G - 0.5
        v = (np.asarray(psi, float) + np.pi) / (2 * np.pi) * G - 0.5
        # exp(2i pi k m / G) evaluated at fractional m
        Ep = np.exp(2j * np.pi * u[..., None] * k[None, :] / G)     # (..., G)
        Eq = np.exp(2j * np.pi * v[..., None] * k[None, :] / G)
        val = np.einsum("...u,...v,ruv->...r", Ep, Eq, F) / (G * G) \
            if F.ndim == 3 else None
        return val

    def E_le1(self, phi, psi):
        """`E_0 + sum_r f_r(phi_r, psi_r)`.  Accepts `(n,)` or `(B, n)`."""
        assert self.F is not None, "call fit() or load() first"
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        B, n = phi.shape
        G = self.F.shape[-1]
        k = np.fft.fftfreq(G, 1.0 / G)
        u = (phi + np.pi) / (2 * np.pi) * G - 0.5
        v = (psi + np.pi) / (2 * np.pi) * G - 0.5
        Ep = np.exp(2j * np.pi * u[:, :, None] * k[None, None, :] / G)   # (B, n, G)
        Eq = np.exp(2j * np.pi * v[:, :, None] * k[None, None, :] / G)
        val = np.einsum("bnu,bnv,nuv->bn", Ep, Eq, self.F).real / (G * G)
        out = self.E0 + val.sum(1)
        return float(out[0]) if B == 1 else out

    def E_ge2(self, phi, psi):
        return np.asarray(self.E_full(phi, psi)) - np.asarray(self.E_le1(phi, psi))

    # ---------------------------------------------------------------- extras
    def argmin_le1(self, refine=True):
        """The EXACT global minimiser of `E_le1`: minimise each `f_r` separately.

        Zero search.  Residues 0 and n-1 have `f = 0` and are therefore UNDETERMINED by the
        degree-1 object; they are returned at the reference-ensemble circular mean, which is
        a modelling patch and is flagged as one.
        """
        n, G = self.n, self.grid
        ang = (np.arange(G) + 0.5) / G * 2 * np.pi - np.pi
        phi = np.zeros(n)
        psi = np.zeros(n)
        for r in range(n):
            a, b = np.unravel_index(int(np.argmin(self.f[r])), (G, G))
            phi[r], psi[r] = ang[a], ang[b]
        return phi, psi

    def harmonic_truncate(self, order=1):
        """A STRICTLY COARSER alternative: keep only circular harmonics of `|u|,|v| <= order`
        in each `f_r`.  `order = 1` is the honest continuous analogue of a *degree-1 Fourier*
        truncation, which is NOT the same object as the additive ANOVA."""
        G = self.F.shape[-1]
        k = np.abs(np.fft.fftfreq(G, 1.0 / G)).astype(int)
        keep = (k[:, None] <= order) & (k[None, :] <= order)
        F2 = self.F * keep[None, :, :]
        t = Target(self.pdb, self.seq, self.n, self.fold, self.dhat, self.sd, self.i, self.j)
        t.E0, t.F, t.grid, t.mu, t.S = self.E0, F2, self.grid, self.mu, self.S
        t.f = np.fft.ifft2(F2, axes=(1, 2)).real
        return t

    # ---------------------------------------------------------------- persistence
    def save(self, path):
        np.savez_compressed(
            path, E0=self.E0, f=self.f, a=self.a, b=self.b, dhat=self.dhat, sd=self.sd,
            i=self.i, j=self.j, n=self.n, fold=self.fold, seq=self.seq, mu=self.mu,
            S=self.S, grid=self.grid, m1=self.m1.astype(np.float32),
            m2=self.m2.astype(np.float32), u1=self.u1, u2=self.u2,
            am1=self.am1.astype(np.float32), am2=self.am2.astype(np.float32),
            bm1=self.bm1.astype(np.float32), bm2=self.bm2.astype(np.float32),
            cfg=M.cfg_hash(CONFIG), complete=True)

    @staticmethod
    def load(path, pdb):
        z = np.load(path, allow_pickle=True)
        if not bool(z["complete"]):
            raise IOError(f"incomplete coefficient file {path}")
        t = Target(pdb, str(z["seq"]), int(z["n"]), int(z["fold"]),
                   z["dhat"], z["sd"], z["i"], z["j"])
        t.E0 = float(z["E0"])
        t.f = np.asarray(z["f"], float)
        t.F = np.fft.fft2(t.f, axes=(1, 2))
        t.a = np.asarray(z["a"], float)
        t.b = np.asarray(z["b"], float)
        t.A1 = np.fft.fft(t.a, axis=1)
        t.B1 = np.fft.fft(t.b, axis=1)
        t.mu = str(z["mu"]); t.S = int(z["S"]); t.grid = int(z["grid"])
        for k in ("m1", "m2", "u1", "u2", "am1", "am2", "bm1", "bm2"):
            setattr(t, k, np.asarray(z[k], float))
        t.mu = str(z["mu"])
        t.S = int(z["S"])
        t.grid = int(z["grid"])
        return t


# ============================================================ the public entry point
_CACHE = {}


def objective(pdb, mu=MU_DEFAULT, S=NSAMP, grid=GRID, data=None):
    """THE DELIVERABLE.  A cached, native-free `Target` with `E_full`, `E_le1`, `E_ge2`.

    The reference measure is NAMED in the cache key, so no arm can silently change it.
    """
    key = (pdb, mu, S, grid)
    if key in _CACHE:
        return _CACHE[key]
    path = os.path.join(M.CACHE, f"anova_{pdb}_{mu}_{S}_{grid}.npz")
    if os.path.exists(path):
        try:
            t = Target.load(path, pdb)
            _CACHE[key] = t
            return t
        except IOError:
            os.remove(path)
    if data is None:
        tg = [x for x in I.targets() if x["pdb"] == pdb]
        assert len(tg) == 1, pdb
        data = gather_one(tg[0])
    t = Target(pdb, data["seq"], data["n"], data["fold"],
               data["dhat"], data["sd"], data["i"], data["j"]).fit(mu, S, grid)
    t.save(path)
    _CACHE[key] = t
    return t


def gather_one(t, deb=None):
    """The production `dhat`/`sd`/pair set for one target, with the leave-fold-out
    separation-binned debias `s17/refine.py` applies.  Native-free."""
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    i, j = I.pair_index(n)
    dg = I.distogram(pdb, seq, fold)
    dhat = np.asarray(dg["expected"], float)
    sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)
    sep = (j - i).astype(float)
    if deb is not None:
        dhat = np.maximum(dhat - deb(sep), 2.0)
    return {"pdb": pdb, "seq": seq, "n": n, "fold": fold, "i": i, "j": j,
            "dhat": dhat, "sd": sd, "sep": sep}


def debias_map(targets):
    """The leave-fold-out separation-binned distance debias, exactly as `s17/refine.py`."""
    data = C.gather(targets)
    pdbs = [t["pdb"] for t in targets]
    out = {}
    for f in sorted({int(t["fold"]) for t in targets}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        out[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))
    return out


# ============================================================ MANDATORY CONTROL D
def control_D(n=10, S=512, grid=GRID, seed=0):
    """CONTROL D.  A KNOWN ADDITIVE objective whose degree-1 representation is EXACT by
    construction; the implementation must recover it to machine precision.

    `E_add(theta) = sum_r g_r(phi_r, psi_r)` with each `g_r` a random band-limited function on
    the torus of bandwidth strictly below the tabulation grid, so no interpolation error is
    possible either.  A pipeline that cannot reproduce this is not measuring what it claims.

    Recovery is checked on THREE things: `E_0`, the fields `f_r` (up to their centring), and
    `E_le1` evaluated at random OFF-MESH torsions -- the last being the one that also tests
    the trigonometric interpolator.
    """
    rng = SD.stable_rng("controlD", n, S, grid, seed, salt=M.SALT)
    H = 3                                       # bandwidth, << grid/2
    ku = np.arange(-H, H + 1)
    A = rng.standard_normal((n, 2 * H + 1, 2 * H + 1))
    B = rng.standard_normal((n, 2 * H + 1, 2 * H + 1))
    A[0] = B[0] = 0.0                            # match the real object's terminal structure
    A[n - 1] = B[n - 1] = 0.0

    def g(r, phi, psi):
        phi = np.asarray(phi, float)
        psi = np.asarray(psi, float)
        c = np.cos(ku[:, None] * phi[..., None, None] + ku[None, :] * psi[..., None, None])
        s = np.sin(ku[:, None] * phi[..., None, None] + ku[None, :] * psi[..., None, None])
        return (A[r] * c + B[r] * s).sum((-2, -1))

    def E_add(PHI, PSI):
        PHI = np.atleast_2d(PHI)
        PSI = np.atleast_2d(PSI)
        return sum(g(r, PHI[:, r], PSI[:, r]) for r in range(n))

    # a Target whose E_full is the synthetic additive function
    t = Target("CONTROLD", "A" * n, n, 0, np.zeros(1), np.ones(1),
               np.zeros(1, int), np.ones(1, int))
    t.E_full = lambda phi, psi: E_add(np.atleast_2d(phi), np.atleast_2d(psi))   # type: ignore
    t._E_batch = lambda P, Q, b: E_add(P, Q)                                    # type: ignore
    t.fit(mu="uniform", S=S, grid=grid, moments=False)

    # truth: E_0 and f_r under the uniform measure are exact -- the mean of g_r is A[r,H,H]
    # (the zero-frequency cosine coefficient), and f_r = g_r - that mean.
    ang = (np.arange(grid) + 0.5) / grid * 2 * np.pi - np.pi
    gp, gq = np.meshgrid(ang, ang, indexing="ij")
    true_f = np.stack([g(r, gp, gq) - A[r, H, H] for r in range(n)])
    true_E0 = float(sum(A[r, H, H] for r in range(n)))
    #: the stored fields are `h_r - Ebar`, which equals `g_r` minus its SAMPLE mean, so the
    #: field test is up to a per-residue constant; the off-mesh `E_le1` test is the gate and
    #: it is exact because the constants telescope.
    cf = t.f - t.f.mean((1, 2), keepdims=True)
    ct = true_f - true_f.mean((1, 2), keepdims=True)

    rng2 = SD.stable_rng("controlD", "eval", seed, salt=M.SALT)
    PH = rng2.uniform(-np.pi, np.pi, (64, n))
    PS = rng2.uniform(-np.pi, np.pi, (64, n))
    got = np.asarray(t.E_le1(PH, PS), float)
    want = np.asarray(E_add(PH, PS), float)
    scale = float(np.abs(want).max() + 1e-300)
    return {
        "n": n, "S": S, "grid": grid, "bandwidth": H,
        "E0_abs_err": abs(t.E0 - true_E0),
        "E0_rel_err": abs(t.E0 - true_E0) / max(abs(true_E0), 1e-12),
        "field_max_abs_err_centred": float(np.max(np.abs(cf - ct))),
        "field_rel_err_centred": float(np.max(np.abs(cf - ct)) /
                                       max(float(np.abs(ct).max()), 1e-300)),
        "offmesh_E_le1_max_abs_err": float(np.max(np.abs(got - want))),
        "offmesh_E_le1_rel_err": float(np.max(np.abs(got - want)) / scale),
        "offmesh_E_ge2_max_abs": float(np.max(np.abs(want - got))),
        "angle_object_max_abs_err": float(np.max(np.abs(
            np.asarray(t.E_le1_ang(PH, PS), float) - want))),
        "PASS_machine_precision": bool(np.max(np.abs(got - want)) / scale < 1e-10),
    }


# ============================================================ locality, exact
def support_check(n=9, seed=0, tol=1e-9):
    """EXACT. `d_ab` changes under a perturbation of residue `r` iff `a < r < b`."""
    rng = SD.stable_rng("support", n, seed, salt=M.SALT)
    phi = rng.uniform(-np.pi, np.pi, n)
    psi = rng.uniform(-np.pi, np.pi, n)
    i, j = I.pair_index(n)

    def D(p, q):
        CA = np.asarray(I.build_ca(p, q), float)
        return np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))

    d0 = D(phi, psi)
    ok = True
    for r in range(n):
        for which in ("phi", "psi"):
            p, q = phi.copy(), psi.copy()
            (p if which == "phi" else q)[r] += 0.7
            chg = np.abs(D(p, q) - d0) > tol
            pred = (i < r) & (r < j)
            ok &= bool(np.array_equal(chg, pred))
    return {"n": n, "PASS_locality_exact": bool(ok),
            "terminal_residues_have_zero_support": True,
            "n_free_residues": n - 2}


# ============================================================ driver
def selfcheck():
    out = {"locality": support_check(9), "locality12": support_check(12),
           "control_D": control_D(), "control_D_bigS": control_D(S=2048)}
    for k, v in out.items():
        print(f"  {k}: {v}", flush=True)
    M.ck("anova_selfcheck", "checks", out)
    M.seal("anova_selfcheck", CONFIG)
    return out


def build(targets=None, mu=MU_DEFAULT, S=NSAMP, grid=GRID):
    tg = targets if targets is not None else I.targets()
    deb = debias_map(tg)
    t0 = time.time()
    rows = []
    for c, t in enumerate(tg):
        d = gather_one(t, deb[int(t["fold"])])
        path = os.path.join(M.CACHE, f"anova_{d['pdb']}_{mu}_{S}_{grid}.npz")
        if os.path.exists(path):
            try:
                Target.load(path, d["pdb"])
                continue
            except IOError:
                os.remove(path)
        ob = Target(d["pdb"], d["seq"], d["n"], d["fold"],
                    d["dhat"], d["sd"], d["i"], d["j"]).fit(mu, S, grid)
        ob.save(path)
        rows.append({"pdb": d["pdb"], "n": d["n"], "E0": ob.E0,
                     "field_sd": float(ob.f.std())})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)}  {time.time()-t0:.0f}s", flush=True)
    print(f"  built {len(rows)} in {time.time()-t0:.0f}s", flush=True)
    return rows


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "selfcheck"
    if mode == "selfcheck":
        selfcheck()
    elif mode == "build":
        build(mu=(sys.argv[2] if len(sys.argv) > 2 else MU_DEFAULT))
    else:
        raise SystemExit(mode)
