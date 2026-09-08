"""SPRINT 15, ALIGN -- shared machinery for ENGINEERING error alignment.

WHY THIS WORKSTREAM EXISTS.

`s15/info_FINDINGS.md` C.3 measured that the i.i.d. phase surface OVER-PRICES every real
torsion channel by 0.6-2.9 A, and named the mechanism: real errors sit in the RMSD-QUIET
directions of torsion space.  The statistic is

    alignment(e) = ||J e|| / (||e|| * s_rms),   s_rms = sqrt(trace(J'J) / 2n)

= damage per unit angular error relative to a random direction (0.945 for a random Gaussian
direction at these dimensions, 0.565 for the best fragment channel).  Every aligned emitter in
that study is aligned because it is a real fragment.  **Nothing measured says a generator can
be steered into the quiet subspace.**  This workstream asks whether it can.

WHAT IS NEW HERE, RELATIVE TO `s15/info_null.py`.

`info_null` computes J by CENTRAL DIFFERENCES (2n forward builds per target) at the NATIVE
torsions.  Both facts matter:

  * The central-difference Jacobian is an O(n) x O(n) cost that is fine for a diagnostic and
    far too slow inside a fit.  This module supplies the **exact analytic** Jacobian from the
    same rigid-suffix-rotation identity that `core.project._torsion_grad` uses in reverse
    mode:  d CA_k / d t_s = u_s x (CA_k - a_s)  for the suffix the torsion moves.
    `align_jac.verify()` checks it against `info_null.jacobian` -- an independent
    implementation, as the brief requires, not a second copy of the same code.

  * Evaluating J at the NATIVE torsions is ORACLE.  Every function here takes the torsions it
    is handed, so the SAME code gives an ORACLE DIAGNOSTIC when fed native torsions and a
    fully NATIVE-FREE quantity when fed the current working structure.  The distinction is
    carried in the caller's arm name, never hidden here.

THE SUPERPOSITION.  RMSD is invariant under rigid motion, so the six rigid-body directions
carry no RMSD and must be projected out of J before any spectrum is read.  At the base point
the Kabsch-superposed derivative is exactly the raw derivative projected orthogonally to the
span of {3 translations, 3 infinitesimal rotations about the centroid}; `align_jac.verify()`
confirms this against the finite-difference-with-Kabsch construction.  With J superposed,

    first-order CA-RMSD produced by an angular error e  =  ||J e|| / sqrt(n).
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from core import project as pj               # noqa: E402
from s12 import instrument as I              # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
CACHE = os.path.join(ROOT, "s15", "cache")
for _d in (RESULTS, CACHE):
    os.makedirs(_d, exist_ok=True)


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


# ------------------------------------------------------------------ the analytic Jacobian
def frames_ca(phi, psi):
    phi = np.asarray(phi, float)
    psi = np.asarray(psi, float)
    G = pj.frames(phi[None], psi[None])[0]
    CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
    return G, CA


def raw_jacobian(phi, psi):
    """`(n, 3, 2n)` d CA / d(phi_0..phi_{n-1}, psi_0..psi_{n-1}), exact, no superposition.

    Column order is exactly `core.project._torsion_grad`'s output order, so a gradient from
    that routine and a column of this Jacobian index the same torsion.

    `phi[0]`, `phi[n-1]` and `psi[n-1]` are never read by the builder, so their columns are
    identically zero here.  `psi[0]` is NOT zero here -- it is a rigid rotation of the whole
    chain about an axis through CA[0] -- and becomes zero only after the superposition
    projection below.  That is the exact content of the "four inert torsions" correction.
    """
    G, CA = frames_ca(phi, psi)
    n = CA.shape[0]
    J = np.zeros((n, 3, 2 * n))
    for j in range(1, n):                     # phi[j] applied at step 3j, frame G[3j-1]
        F = G[3 * j - 1]
        J[j + 1:, :, j] = np.cross(F[:3, 0], CA[j + 1:] - F[:3, 3])
    for i in range(0, n - 1):                 # psi[i] applied at step 3i+1, frame G[3i]
        F = G[3 * i]
        J[i + 1:, :, n + i] = np.cross(F[:3, 0], CA[i + 1:] - F[:3, 3])
    return J, CA


def rigid_basis(CA):
    """Orthonormal `(3n, 6)` basis of the rigid-body tangent directions at `CA`."""
    n = CA.shape[0]
    c = CA.mean(0)
    cols = []
    for a in range(3):
        m = np.zeros((n, 3))
        m[:, a] = 1.0
        cols.append(m.ravel())
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        cols.append(np.cross(e, CA - c).ravel())
    B = np.column_stack(cols)
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    return U[:, s > 1e-9 * s[0]]


def sup_jacobian(phi, psi):
    """`(J_sup (3n,2n), J_raw (n,3,2n), CA (n,3))`.  `||J_sup e|| / sqrt(n)` is the RMSD."""
    Jr, CA = raw_jacobian(phi, psi)
    n = CA.shape[0]
    J = Jr.reshape(3 * n, 2 * n)
    Q = rigid_basis(CA)
    return J - Q @ (Q.T @ J), Jr, CA


def s_rms(J):
    """Root-mean-square damage of a random unit torsion direction through `J`."""
    m = J.shape[1]
    return float(np.sqrt((J ** 2).sum() / m))


def alignment(J, e):
    """`||J e|| / (||e|| * s_rms)` -- the INFO workstream's statistic, same definition."""
    e = np.asarray(e, float)
    ne = float(np.linalg.norm(e))
    if ne < 1e-12:
        return float("nan")
    return float(np.linalg.norm(J @ e) / (ne * s_rms(J)))


# ------------------------------------------------------------------ restraint sensitivities
def pair_jacobian(Jr, CA, i, j):
    """`(P, 2n)` d d_ij / d(theta) for the restraint pairs `(i, j)`.

    Distances are rigid-invariant, so the RAW Jacobian is used; projecting first changes
    nothing and costs a matmul.
    """
    r = CA[i] - CA[j]
    d = np.maximum(np.sqrt((r * r).sum(1)), 1e-9)
    u = r / d[:, None]
    return np.einsum("pa,pak->pk", u, Jr[i] - Jr[j]), d


def pair_response(phi, psi, i, j):
    """`(P,)` ||d d_ij / d theta|| -- how strongly the structure responds to each restraint.

    NATIVE-FREE when evaluated at a working structure.  This is the quantity the
    Jacobian-weighted arm reweights by.
    """
    Jr, CA = raw_jacobian(phi, psi)
    G, _ = pair_jacobian(Jr, CA, i, j)
    return np.linalg.norm(G, axis=1)


# ------------------------------------------------------------------ subspace bookkeeping
def spectrum(J):
    """`(s (2n,), V (2n, 2n))` singular values (descending) and right singular vectors."""
    _, s, Vt = np.linalg.svd(J, full_matrices=True)
    m = J.shape[1]
    s = np.concatenate([s, np.zeros(m - len(s))]) if len(s) < m else s
    return s, Vt.T


def energy_split(V, s, e, thresh):
    """Fraction of ||e||^2 and of ||J e||^2 carried by directions with `s_i <= thresh`."""
    c = V.T @ np.asarray(e, float)
    q = c * c
    dmg = (s * s) * q
    lo = s <= thresh
    tq, td = q.sum(), dmg.sum()
    return (float(q[lo].sum() / tq) if tq > 0 else float("nan"),
            float(dmg[lo].sum() / td) if td > 0 else float("nan"))


# ------------------------------------------------------------------ the fit
def fit(dhat, sd, i, j, phi0, psi0, wpair=None, kind="squared", s=1.0,
        basis=None, theta0=None, maxiter=400):
    """One restraint fit.  A superset of `s15.robust.fit_robust`.

    `wpair`   extra multiplicative per-restraint weight on top of 1/sd^2 (the alignment arms)
    `kind`/`s`  the robust loss and its scale, identical definitions to `s15.robust.rho`
    `basis`   `(2n, r)` -- if given, optimise over `theta = theta0 + basis @ c` only, which
              is the "projection onto the well-conditioned subspace" arm.  The chain rule is
              one extra matmul, so the gradient stays exact.

    Returns `(phi, psi, f)`.  Selection among starts is by `f` -- the OBJECTIVE -- only.
    """
    n = len(phi0)
    inv = 1.0 / np.asarray(sd, float)
    w = np.ones(len(dhat)) if wpair is None else np.asarray(wpair, float)
    from s15.robust import rho

    def fg_theta(x):
        phi = x[:n]
        psi = x[n:]
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[i] - CA[j]
        d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        r = (d - dhat) * inv
        v, g = rho(r, kind, s)
        f = float((w * v).sum())
        coef = ((w * g * inv) / d)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, i, coef)
        np.add.at(gCA, j, -coef)
        return f, pj._torsion_grad(G, CA, gCA)

    from scipy.optimize import minimize
    if basis is None:
        x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
        res = minimize(fg_theta, x0, jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12,
                                "gtol": 1e-10})
        return res.x[:n], res.x[n:], float(res.fun)

    B = np.asarray(basis, float)
    t0 = (np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
          if theta0 is None else np.asarray(theta0, float))

    def fg_c(c):
        f, g = fg_theta(t0 + B @ c)
        return f, B.T @ g

    res = minimize(fg_c, np.zeros(B.shape[1]), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    x = t0 + B @ res.x
    return x[:n], x[n:], float(res.fun)


# ------------------------------------------------------------------ post-hoc scoring
def score_axes(phi, psi, nphi, npsi, nat, dhat, sd, i, j, dtrue):
    """BOTH axes for one emitted structure, so no arm can be read on one only.

    RAW axes (must move if the intervention improved the channel's accuracy):
      `tors_rms_deg`   circular RMS angular error against the native torsions   [ORACLE eval]
      `dist_mae`       MAE of the emitted CA-CA distances against the true ones [ORACLE eval]
      `resid_medz`     median |d - dhat| / sd, the NATIVE-FREE restraint residual
    ALIGNMENT axis:
      `align`          ||J e|| / (||e|| s_rms) with J at the NATIVE torsions    [ORACLE]
      `align_self`     the same with J at the EMITTED torsions -- NATIVE-FREE J, ORACLE e
    OUTCOME:
      `rmsd`           full-chain CA-RMSD, the frozen definition
    """
    ca = np.asarray(I.build_ca(phi, psi), float)
    d = np.sqrt(((ca[i] - ca[j]) ** 2).sum(1))
    e = np.concatenate([wrap(np.asarray(phi, float) - nphi),
                        wrap(np.asarray(psi, float) - npsi)])
    Jn, _, _ = sup_jacobian(nphi, npsi)
    Js, _, _ = sup_jacobian(phi, psi)
    return {"rmsd": float(I.ca_rmsd(ca, nat)),
            "tors_rms_deg": float(np.degrees(np.sqrt((e ** 2).mean()))),
            "dist_mae": float(np.abs(d - dtrue).mean()),
            "resid_medz": float(np.median(np.abs(d - dhat) / np.asarray(sd, float))),
            "align": alignment(Jn, e),
            "align_self": alignment(Js, e)}
