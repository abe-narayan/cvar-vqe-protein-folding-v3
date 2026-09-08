"""SPRINT 18 / MATH -- the adapter the EXPERIMENT workstream consumes.

`s18/COORD_exp_to_math.md` publishes the interface EXPERIMENT will call.  This module
implements it on top of `s18/math_anova.py` so there is ONE definition of the degree-1 object
in the sprint.

------------------------------------------------------------------------------------------
A NAMING CORRECTION THAT IS NOT COSMETIC -- READ THIS BEFORE WIRING ANYTHING UP
------------------------------------------------------------------------------------------

EXPERIMENT's provisional file calls the ANGLE-ADDITIVE object `E_le1` ("STRICT weight-<=1")
and the RESIDUE-ADDITIVE object `E_res`.  The lattice measurement says that split is where
the whole sprint turns, and that the names invite the wrong reading:

  * The BRIEF's bridge -- first-order functional ANOVA under `mu = tensor_i mu_i` -- is an
    ANOVA over whatever product factorisation you name.  With the factor = RESIDUE it is the
    residue-additive object.  That is the object the bridge derives, and on the enumerated
    lattice it is worth **-0.004 A [-0.380, +0.307], 3W/5L/11 ties** against the full
    objective.  It is NOT the 2.411 A result.
  * The 2.411 A result is the STRICT WALSH WEIGHT-<=1 projection, i.e. the ANOVA under the
    QUBIT factorisation.  The two lattice qubits of a residue index k-means clusters of the
    joint (phi, psi) library; they are not phi and psi.  So the angle-additive continuous
    object is **not** the continuous image of the strict lattice object.  There is no such
    image: the lattice has no phi/psi factorisation and the continuum has no qubits.
  * And the strict lattice object is not invariant to RELABELLING THE FOUR TORSION STATES --
    a pure bookkeeping change.  Its argmin over 24 random relabellings averages 2.862 +-
    0.135 A, worse than the full objective's 2.661, with the shipped labelling at the 0th
    percentile.

So both objects should still be run -- the angle-additive one is a legitimate, strictly
coarser truncation worth measuring in its own right -- but it must be labelled as an ANALOGY,
never as the continuous form of the 2.411 A result, and `E_res` is the object the BRIEF's
derivation actually specifies.

To avoid a silent mismatch this adapter exposes BOTH under EXPERIMENT's names AND under
unambiguous ones:

    obj.E_res / obj.E_le1_residue   residue-additive   <-- the BRIEF's ANOVA object
    obj.E_le1 / obj.E_le1_angle     angle-additive     <-- strictly coarser, an analogy

Every method returns `(value, gradient_2n)` with the gradient in the order
`[dphi_0..dphi_{n-1}, dpsi_0..dpsi_{n-1}]`, matching `core.project._torsion_grad`, so
L-BFGS treats every arm identically.  `verify_gradients` checks all three against central
differences.

FOUR STRUCTURAL FACTS EXPERIMENT SHOULD BUILD INTO THE HARNESS (all EXACT):

 1. The distance objective is invariant to `phi_0, psi_0, phi_{n-1}, psi_{n-1}`, so both
    truncations have IDENTICALLY ZERO field on the terminal residues.  `argmin_le1()` cannot
    place them, and the frozen metric scores them.  Whatever they are set to is a modelling
    patch and must be reported as one; the adapter leaves them at the start value.
 2. Both truncations are separable, so their global argmin is closed form.  No search.
 3. `E_lambda = (1-lam) E_le1 + lam E_full`, so `lam = 1` reproduces the full objective
    EXACTLY -- `check_lambda_identity` asserts it to machine precision.
 4. Because `E = sum_p w_p (d_p^2 - 2 dhat_p d_p + dhat_p^2)` is affine in `(d, d^2)` at
    fixed theta, the truncation's tables are affine in `dhat`.  `with_dhat` therefore
    rebuilds either truncated objective for an arbitrary distance vector with NO new chain
    builds -- which is what makes the alpha-ladder and the argmin-sensitivity analysis cheap.
"""
from __future__ import annotations

import numpy as np

from s18 import math_anova as A
from s18 import math_lib as M

from core import project as pj                       # noqa: E402
from s12 import instrument as I                      # noqa: E402


class Obj:
    """EXPERIMENT's interface over one `math_anova.Target`."""

    def __init__(self, t):
        self.t = t
        self.n = t.n
        self.E0 = t.E0
        self.pdb = t.pdb
        self.mu = t.mu
        self.grid = t.grid
        self.S = t.S

    # ------------------------------------------------------------------ full
    def E_full(self, phi, psi):
        """`(f, grad)` for the deployed objective, gradient by the exact torsion chain rule
        (`core.project._torsion_grad`) -- the same path `s15.align_lib.fit` uses, so an arm
        optimised here is optimised identically to `s17/refine.py`."""
        t = self.t
        phi = np.asarray(phi, float).ravel()
        psi = np.asarray(psi, float).ravel()
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[t.i] - CA[t.j]
        d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        r = d - t.dhat
        f = float((t.w * r * r).sum())
        coef = ((2.0 * t.w * r) / d)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, t.i, coef)
        np.add.at(gCA, t.j, -coef)
        return f, np.asarray(pj._torsion_grad(G, CA, gCA), float)

    # ------------------------------------------------------------------ residue-additive
    def E_res(self, phi, psi):
        """THE BRIEF'S OBJECT: `E_0 + sum_r f_r(phi_r, psi_r)`, exact trig interpolation."""
        t = self.t
        phi = np.asarray(phi, float).ravel()
        psi = np.asarray(psi, float).ravel()
        G = t.grid
        k = np.fft.fftfreq(G, 1.0 / G)
        u = (phi + np.pi) / (2 * np.pi) * G - 0.5
        v = (psi + np.pi) / (2 * np.pi) * G - 0.5
        Ep = np.exp(2j * np.pi * u[:, None] * k[None, :] / G)          # (n, G)
        Eq = np.exp(2j * np.pi * v[:, None] * k[None, :] / G)
        c = 2j * np.pi / (2 * np.pi) * 1.0                              # d/dangle of exp
        dEp = Ep * (1j * k[None, :])
        dEq = Eq * (1j * k[None, :])
        val = np.einsum("nu,nv,nuv->n", Ep, Eq, t.F).real / (G * G)
        gp = np.einsum("nu,nv,nuv->n", dEp, Eq, t.F).real / (G * G)
        gq = np.einsum("nu,nv,nuv->n", Ep, dEq, t.F).real / (G * G)
        del c
        return float(t.E0 + val.sum()), np.concatenate([gp, gq])

    E_le1_residue = E_res

    # ------------------------------------------------------------------ angle-additive
    def E_ang(self, phi, psi):
        """The strictly coarser sub-residue truncation:
        `E_0 + sum_r a_r(phi_r) + sum_r b_r(psi_r)`."""
        t = self.t
        phi = np.asarray(phi, float).ravel()
        psi = np.asarray(psi, float).ravel()
        G = t.grid
        k = np.fft.fftfreq(G, 1.0 / G)
        u = (phi + np.pi) / (2 * np.pi) * G - 0.5
        v = (psi + np.pi) / (2 * np.pi) * G - 0.5
        Ep = np.exp(2j * np.pi * u[:, None] * k[None, :] / G)
        Eq = np.exp(2j * np.pi * v[:, None] * k[None, :] / G)
        va = np.einsum("nu,nu->n", Ep, t.A1).real / G
        vb = np.einsum("nu,nu->n", Eq, t.B1).real / G
        ga = np.einsum("nu,nu->n", Ep * (1j * k[None, :]), t.A1).real / G
        gb = np.einsum("nu,nu->n", Eq * (1j * k[None, :]), t.B1).real / G
        return float(t.E0 + va.sum() + vb.sum()), np.concatenate([ga, gb])

    E_le1 = E_ang
    E_le1_angle = E_ang

    # ------------------------------------------------------------------ complements
    def E_ge2(self, phi, psi, base="res"):
        """`E_full - E_trunc`, with its gradient.  `base` names WHICH truncation."""
        ff, gf = self.E_full(phi, psi)
        lf, gl = (self.E_res if base == "res" else self.E_ang)(phi, psi)
        return ff - lf, gf - gl

    def E_lambda(self, phi, psi, lam, base="res"):
        """`(1-lam) E_trunc + lam E_full`.  `lam = 1` is the full objective EXACTLY."""
        ff, gf = self.E_full(phi, psi)
        lf, gl = (self.E_res if base == "res" else self.E_ang)(phi, psi)
        return (1 - lam) * lf + lam * ff, (1 - lam) * gl + lam * gf

    # ------------------------------------------------------------------ exact argmins
    def argmin_res(self, phi0=None, psi0=None):
        """EXACT global minimiser of the residue-additive object -- coordinate-wise, no
        search.  Terminal residues are undetermined (zero field) and are LEFT AT THE START;
        pass `phi0/psi0` or they default to the reference-ensemble mode."""
        t = self.t
        G = t.grid
        ang = (np.arange(G) + 0.5) / G * 2 * np.pi - np.pi
        phi = np.zeros(t.n) if phi0 is None else np.asarray(phi0, float).copy()
        psi = np.zeros(t.n) if psi0 is None else np.asarray(psi0, float).copy()
        for r in range(1, t.n - 1):
            a, b = np.unravel_index(int(np.argmin(t.f[r])), (G, G))
            phi[r], psi[r] = ang[a], ang[b]
        return phi, psi

    def argmin_ang(self, phi0=None, psi0=None):
        t = self.t
        G = t.grid
        ang = (np.arange(G) + 0.5) / G * 2 * np.pi - np.pi
        phi = np.zeros(t.n) if phi0 is None else np.asarray(phi0, float).copy()
        psi = np.zeros(t.n) if psi0 is None else np.asarray(psi0, float).copy()
        for r in range(1, t.n - 1):
            phi[r] = ang[int(np.argmin(t.a[r]))]
            psi[r] = ang[int(np.argmin(t.b[r]))]
        return phi, psi

    argmin_le1 = argmin_ang

    # ------------------------------------------------------------------ dhat handle
    def with_dhat(self, dhat):
        """A new `Obj` for a different distance vector.  EXACT, no new chain builds."""
        return Obj(self.t.with_dhat(dhat))

    @property
    def terminal_residues_undetermined(self):
        return (0, self.n - 1)


def build(pdb=None, mu=A.MU_DEFAULT, S=A.NSAMP, grid=A.GRID, target=None):
    """THE ENTRY POINT.  `build('1A13')` returns a cached, native-free `Obj`.

    `target` may be a `s12.instrument.targets()` row if the caller already has it.
    """
    if target is not None and pdb is None:
        pdb = target["pdb"]
    t = A.objective(pdb, mu=mu, S=S, grid=grid,
                    data=A.gather_one(target) if target is not None else None)
    return Obj(t)


# ====================================================================== verification
def verify_gradients(pdb, eps=1e-6, seed=0):
    """Central-difference check of all three gradients.  Run before trusting any optimiser."""
    o = build(pdb)
    from s15 import seed as SD
    rng = SD.stable_rng(pdb, seed, "gradcheck", salt=M.SALT)
    phi = rng.uniform(-np.pi, np.pi, o.n)
    psi = rng.uniform(-np.pi, np.pi, o.n)
    out = {}
    for name, fn in (("E_full", o.E_full), ("E_res", o.E_res), ("E_ang", o.E_ang)):
        f0, g = fn(phi, psi)
        num = np.empty(2 * o.n)
        for c in range(2 * o.n):
            for s in (+1, -1):
                p, q = phi.copy(), psi.copy()
                (p if c < o.n else q)[c % o.n] += s * eps
                if s > 0:
                    fp = fn(p, q)[0]
                else:
                    fm = fn(p, q)[0]
            num[c] = (fp - fm) / (2 * eps)
        den = max(float(np.abs(num).max()), 1e-12)
        out[name] = {"max_abs_err": float(np.abs(num - g).max()),
                     "rel_err": float(np.abs(num - g).max() / den)}
    return out


def check_lambda_identity(pdb, seed=0):
    """EXACT. `E_lambda(., 1) == E_full` and `E_lambda(., 0) == E_trunc`."""
    o = build(pdb)
    from s15 import seed as SD
    rng = SD.stable_rng(pdb, seed, "lamcheck", salt=M.SALT)
    phi = rng.uniform(-np.pi, np.pi, o.n)
    psi = rng.uniform(-np.pi, np.pi, o.n)
    f1, _ = o.E_lambda(phi, psi, 1.0)
    f0, _ = o.E_lambda(phi, psi, 0.0)
    return {"lam1_minus_full": float(f1 - o.E_full(phi, psi)[0]),
            "lam0_minus_res": float(f0 - o.E_res(phi, psi)[0])}


def check_refine_equivalence(pdb, seed=0):
    """`E_full` here must equal `s17/refine.py::_obj` at the same torsions, to the last bit."""
    o = build(pdb)
    from s15 import seed as SD
    rng = SD.stable_rng(pdb, seed, "refeq", salt=M.SALT)
    phi = rng.uniform(-np.pi, np.pi, o.n)
    psi = rng.uniform(-np.pi, np.pi, o.n)
    CA = I.build_ca(phi, psi)
    d = np.sqrt(((CA[o.t.i] - CA[o.t.j]) ** 2).sum(1))
    ref = float((((d - o.t.dhat) / o.t.sd) ** 2).sum())
    return {"E_full_minus_s17_obj": float(o.E_full(phi, psi)[0] - ref)}


if __name__ == "__main__":
    import sys
    pdb = sys.argv[1] if len(sys.argv) > 1 else "1A13"
    print("gradients        ", verify_gradients(pdb))
    print("lambda identity  ", check_lambda_identity(pdb))
    print("s17 objective eq ", check_refine_equivalence(pdb))
