"""s18/exp_obj.py -- the EXPERIMENT workstream's adapter onto MATH's deliverable.

THE RULE THIS FILE EXISTS TO OBEY.  The BRIEF forbids a second definition of the degree-1
object.  So the objective, the reference measure, the tabulation grid and the quadrature all
come from `s18/math_anova.py` and nothing here redefines them.  What this file adds is the two
things MATH's `Target` does not carry and the experiment cannot run without:

  1. ANALYTIC GRADIENTS of `E_full` and `E_le1`, so every arm is optimised by the same
     L-BFGS to the same tolerance.  A comparison in which one objective is optimised better
     than the other measures the optimiser, not the objective.  Both gradients are verified
     against central differences OF MATH'S OWN CALLABLES in `verify()` -- not against a
     re-derivation, which would only prove this file is self-consistent.
  2. The lambda ladder `E_lam = E_le1 + lam (E_full - E_le1) = (1-lam) E_le1 + lam E_full`,
     which is a convex mix of MATH's two functions and introduces nothing new.

Also here, because they are EXPERIMENT decisions rather than MATH ones:

  * `argmin_le1_hold` -- MATH proves `f_0 = f_{n-1} = 0` EXACTLY, so the degree-1 object
    leaves the two terminal residues completely undetermined while the frozen metric scores
    them.  MATH's `argmin_le1` resolves that by taking `argmin` of an all-zero table, which
    silently returns mesh cell (0, 0), i.e. a fixed corner of the torus.  That is an arbitrary
    structural commitment on 2 of n residues and it would be charged to the degree-1 object as
    if it were information.  This adapter keeps BOTH: MATH's as shipped, and a variant that
    holds the undetermined residues at the START value.  The difference between them prices
    exactly how much of any degree-1 result is the terminal-residue patch.
  * `rand_field` -- the zero-information control an ADDITIVE object needs: fields with the
    same per-residue amplitude and the same band-limit, carrying no information at all.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import scipy.optimize as SO                    # noqa: E402

from core import project as pj                 # noqa: E402
from s18 import math_anova as MA               # noqa: E402


class Obj:
    """Gradient-carrying wrapper around one `math_anova.Target`.  Adds no new definition."""

    def __init__(self, t):
        self.t = t
        self.n = t.n
        self.i, self.j = t.i, t.j
        self.dhat, self.w = t.dhat, t.w
        self.G = t.F.shape[-1]
        self.k = np.fft.fftfreq(self.G, 1.0 / self.G)
        self.nev = {"full": 0, "le1": 0}

    # ------------------------------------------------------------------ full objective
    def E_full(self, phi, psi):
        phi = np.asarray(phi, float)
        psi = np.asarray(psi, float)
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[self.i] - CA[self.j]
        d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        r = d - self.dhat
        f = float((self.w * r * r).sum())
        coef = ((2.0 * self.w * r) / d)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, self.i, coef)
        np.add.at(gCA, self.j, -coef)
        self.nev["full"] += 1
        return f, pj._torsion_grad(G, CA, gCA)

    # ------------------------------------------------------------------ degree-1 object
    def E_le1(self, phi, psi, F=None, E0=None):
        """`E_0 + sum_r f_r(phi_r, psi_r)` and its analytic gradient.

        Same interpolant as `math_anova.Target.E_le1`: the mesh point of angle `x` is
        `m = (x + pi)/(2 pi) G - 0.5`, and `d/dx exp(2 pi i u m / G) = i u exp(...)`, so the
        gradient is one extra einsum over the same tensors.
        """
        F = self.t.F if F is None else F
        E0 = self.t.E0 if E0 is None else E0
        G, k = self.G, self.k
        u = (np.asarray(phi, float) + np.pi) / (2 * np.pi) * G - 0.5
        v = (np.asarray(psi, float) + np.pi) / (2 * np.pi) * G - 0.5
        Ep = np.exp(2j * np.pi * u[:, None] * k[None, :] / G)      # (n, G)
        Eq = np.exp(2j * np.pi * v[:, None] * k[None, :] / G)
        val = np.einsum("nu,nv,nuv->n", Ep, Eq, F) / (G * G)
        gp = np.einsum("nu,u,nv,nuv->n", Ep, 1j * k, Eq, F).real / (G * G)
        gs = np.einsum("nu,nv,v,nuv->n", Ep, Eq, 1j * k, F).real / (G * G)
        self.nev["le1"] += 1
        return float(E0 + val.real.sum()), np.concatenate([gp, gs])

    def E_lambda(self, phi, psi, lam):
        if lam == 0.0:
            return self.E_le1(phi, psi)
        if lam == 1.0:
            return self.E_full(phi, psi)
        a, ga = self.E_le1(phi, psi)
        b, gb = self.E_full(phi, psi)
        return (1 - lam) * a + lam * b, (1 - lam) * ga + lam * gb

    # ------------------------------------------------------------------ certified optima
    def argmin_le1(self, F=None, hold=None):
        """The EXACT global minimiser of the additive object.  Separability, not search.

        `hold` (a `(phi0, psi0)` pair) is used for residues whose field is identically zero --
        by MATH's locality theorem that is residues 0 and n-1 on every target.  With
        `hold=None` this reproduces MATH's shipped behaviour (mesh cell 0, a torus corner).
        """
        F = self.t.F if F is None else F
        f = np.fft.ifft2(F, axes=(1, 2)).real
        G = self.G
        ang = (np.arange(G) + 0.5) / G * 2 * np.pi - np.pi
        ph = np.zeros(self.n)
        ps = np.zeros(self.n)
        for r in range(self.n):
            if hold is not None and np.ptp(f[r]) < 1e-12:
                ph[r], ps[r] = hold[0][r], hold[1][r]
                continue
            a, b = np.unravel_index(int(np.argmin(f[r])), (G, G))
            ph[r], ps[r] = ang[a], ang[b]
        return ph, ps

    def undetermined(self):
        """Which residues the degree-1 object exerts no force on at all."""
        return np.where(np.ptp(self.t.f.reshape(self.n, -1), axis=1) < 1e-12)[0]

    # ------------------------------------------------------------------ the controls
    def rand_field(self, rng):
        """Random additive fields with MATCHED per-residue amplitude and band-limit.

        The zero-information reference for an ADDITIVE object: same functional form, same
        smoothness, same magnitude, no information about the target.  If refining toward THIS
        moves RMSD as much as refining toward `E_le1`, degree-1 carries nothing.
        """
        f = self.t.f
        G = self.G
        amp = f.reshape(self.n, -1).std(1)
        R = rng.standard_normal((self.n, G, G))
        FR = np.fft.fft2(R, axes=(1, 2))
        #: match the radial power spectrum of the real fields, so smoothness matches too
        P = np.abs(np.fft.fft2(f, axes=(1, 2)))
        FR = FR / np.maximum(np.abs(FR), 1e-12) * P
        r2 = np.fft.ifft2(FR, axes=(1, 2)).real
        r2 -= r2.reshape(self.n, -1).mean(1)[:, None, None]
        s = r2.reshape(self.n, -1).std(1)
        r2 *= (amp / np.maximum(s, 1e-12))[:, None, None]
        r2[amp < 1e-12] = 0.0
        return np.fft.fft2(r2, axes=(1, 2))


class ObjAng:
    """Gradient wrapper on MATH's SUB-RESIDUE object `E_0 + sum a_r(phi_r) + sum b_r(psi_r)`.

    Strictly coarser than the residue-additive object, and the closest legitimate continuous
    relative of the lattice's per-qubit truncation.  MATH's own finding is that NO continuous
    object equals the strict Walsh weight-<=1 projection -- the two lattice bits index k-means
    clusters of the joint (phi, psi) library, not phi and psi -- which is the content of F5.
    This arm is therefore reported as "the closest continuous relative", never as "the same
    object", and it exists so that the per-qubit / per-residue distinction is MEASURED rather
    than argued.
    """

    def __init__(self, t):
        self.t = t
        self.n = t.n
        self.G = t.A1.shape[-1]
        self.k = np.fft.fftfreq(self.G, 1.0 / self.G)
        self.nev = 0

    def _one(self, Fk, ang):
        G, k = self.G, self.k
        u = (np.asarray(ang, float) + np.pi) / (2 * np.pi) * G - 0.5
        E = np.exp(2j * np.pi * u[:, None] * k[None, :] / G)
        v = np.einsum("nu,nu->n", E, Fk).real / G
        d = np.einsum("nu,u,nu->n", E, 1j * k, Fk).real / G
        return v, d

    def E(self, phi, psi):
        va, ga = self._one(self.t.A1, np.asarray(phi, float))
        vb, gb = self._one(self.t.B1, np.asarray(psi, float))
        self.nev += 1
        return float(self.t.E0 + va.sum() + vb.sum()), np.concatenate([ga, gb])

    def argmin(self, hold=None):
        """EXACT: each angle independently.  Undetermined angles are held at the start."""
        G = self.G
        ang = (np.arange(G) + 0.5) / G * 2 * np.pi - np.pi
        out = []
        for tab, h in ((self.t.a, None if hold is None else hold[0]),
                       (self.t.b, None if hold is None else hold[1])):
            x = np.empty(self.n)
            for r in range(self.n):
                if h is not None and np.ptp(tab[r]) < 1e-12:
                    x[r] = h[r]
                else:
                    x[r] = ang[int(np.argmin(tab[r]))]
            out.append(x)
        return out[0], out[1]


# ------------------------------------------------------------------------------ optimisation
def lbfgs(fun, phi0, psi0, maxiter=400):
    n = len(phi0)

    def fg(x):
        f, g = fun(x[:n], x[n:])
        return f, np.asarray(g, float)

    r = SO.minimize(fg, np.concatenate([phi0, psi0]), jac=True, method="L-BFGS-B",
                    options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return r.x[:n], r.x[n:], float(r.fun), int(r.nfev), int(r.nit)


# ------------------------------------------------------------------------------ verification
def verify(pdb="1A13", S=256, grid=24):
    """Gradients against central differences OF MATH'S OWN CALLABLES, plus the identity
    `E_le1 + E_ge2 == E_full` read through MATH's functions and mine."""
    from s12 import instrument as I
    from s15 import seed as SD
    tg = I.targets()
    t = [x for x in tg if x["pdb"] == pdb][0]
    d = MA.gather_one(t, MA.debias_map(tg)[int(t["fold"])])
    T = MA.Target(d["pdb"], d["seq"], d["n"], d["fold"],
                  d["dhat"], d["sd"], d["i"], d["j"]).fit("pool", S, grid)
    o = Obj(T)
    rng = SD.stable_rng("verify", pdb)
    n = T.n
    phi = rng.uniform(-np.pi, np.pi, n)
    psi = rng.uniform(-np.pi, np.pi, n)
    th = np.concatenate([phi, psi])
    out = {}
    for nm, mine, theirs in (("full", o.E_full, T.E_full), ("le1", o.E_le1, T.E_le1)):
        f, g = mine(phi, psi)
        out[f"{nm}_value_vs_math"] = float(abs(f - float(theirs(phi, psi))))
        fd = np.empty(2 * n)
        for k in range(2 * n):
            e = np.zeros(2 * n); e[k] = 1e-5
            fd[k] = (float(theirs((th + e)[:n], (th + e)[n:]))
                     - float(theirs((th - e)[:n], (th - e)[n:]))) / 2e-5
        out[f"{nm}_grad_maxerr_vs_math_FD"] = float(np.abs(g - fd).max())
        out[f"{nm}_grad_relerr"] = float(np.abs(g - fd).max()
                                         / max(np.abs(fd).max(), 1e-12))
    out["le1_plus_ge2_minus_full"] = float(abs(
        o.E_le1(phi, psi)[0] + float(T.E_ge2(phi, psi)) - o.E_full(phi, psi)[0]))
    out["undetermined_residues"] = o.undetermined().tolist()
    ap, aq = o.argmin_le1()
    bp, bq = o.argmin_le1(hold=(phi, psi))
    out["argmin_math_le1"] = float(T.E_le1(ap, aq))
    out["argmin_hold_le1"] = float(T.E_le1(bp, bq))
    out["argmin_identical_off_terminals"] = bool(
        np.allclose(ap[1:n - 1], bp[1:n - 1]) and np.allclose(aq[1:n - 1], bq[1:n - 1]))
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(verify(), indent=1))
