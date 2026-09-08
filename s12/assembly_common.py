"""Sprint 12 / assembly: shared primitives (piece placement, torsion concatenation, junction
optimisation, cut enumeration).  Nothing here reads a native except functions whose names
start with `oracle_`.
"""
from __future__ import annotations
import os, sys, math, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import assembly_bank as AB

LMIN, LMAX = 3, 10


def cuts_2(n):
    """All 2-piece cut points c (A=[0,c), B=[c,n)) with both lengths in the bank range."""
    return [c for c in range(LMIN, n - LMIN + 1) if c <= LMAX and n - c <= LMAX]


def cuts_m(n, m, lmin=LMIN, lmax=LMAX):
    """All m-piece length compositions (L1..Lm) of n with lmin <= Li <= lmax."""
    out = []
    for comp in itertools.product(range(lmin, lmax + 1), repeat=m):
        if sum(comp) == n:
            out.append(tuple(comp))
    return out


def equal_split(n, m):
    base, r = divmod(n, m)
    return tuple(base + (1 if k < r else 0) for k in range(m))


def intervals(comp):
    s = 0; out = []
    for L in comp:
        out.append((s, L)); s += L
    return out


def circ_mean(a, b):
    return np.arctan2((np.sin(a) + np.sin(b)) / 2.0, (np.cos(a) + np.cos(b)) / 2.0)


def concat_torsions(n, placed, overlap=0, conv="cut"):
    """placed: list of (s, L, phi(L,), psi(L,)) covering [0,n) as consecutive intervals.
    overlap=1 with conv in {"A","B","avg"}: each piece except the first is extended one
    residue to the LEFT (its phi/psi arrays are length L+1 starting at s-1); the shared
    residue takes torsions from the left piece (A), right piece (B), or a circular mean."""
    phi = np.zeros(n); psi = np.zeros(n)
    for k, (s, L, ph, ps) in enumerate(placed):
        if overlap == 0 or k == 0:
            phi[s:s + L] = ph; psi[s:s + L] = ps
        else:
            # ph/ps are length L+1 starting at s-1
            phi[s:s + L] = ph[1:]; psi[s:s + L] = ps[1:]
            if conv == "B":
                phi[s - 1] = ph[0]; psi[s - 1] = ps[0]
            elif conv == "avg":
                phi[s - 1] = circ_mean(phi[s - 1], ph[0]); psi[s - 1] = circ_mean(psi[s - 1], ps[0])
    return phi, psi


def build_many(PHI, PSI, chunk=20000):
    PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    out = np.zeros((len(PHI), PHI.shape[1], 3))
    for a in range(0, len(PHI), chunk):
        out[a:a + chunk] = I.build_ca(PHI[a:a + chunk], PSI[a:a + chunk])
    return out


def local_rmsd_all(fold, s, L, T, real=False):
    """CA-RMSD of every bank piece (fold, L) to the (L,3) trace T (ideal-built pieces unless real)."""
    P = AB.pieces(fold, L)
    return I.kabsch_rmsd_batch(P["rca" if real else "ca"], T)


def junction_params(comp):
    """Indices of the CA-relevant junction torsions for a composition: (psi index, phi index)."""
    out = []; s = 0
    for L in comp[:-1]:
        s += L; out.append((s - 1, s))
    return out


def optimise_junctions(phi, psi, comp, target, grid=24, refine=True):
    """ORACLE placement floor: optimise only (psi_b, phi_{b+1}) at each junction to minimise the
    global CA-RMSD to `target`.  Coordinate-wise grid then Powell refinement.  Returns
    (best_rmsd, phi, psi)."""
    from scipy.optimize import minimize
    phi = phi.copy(); psi = psi.copy(); J = junction_params(comp)
    ang = np.linspace(-math.pi, math.pi, grid, endpoint=False)

    def rm(ph, ps):
        return I.ca_rmsd(I.build_ca(ph[None], ps[None])[0], target)

    best = rm(phi, psi)
    for _ in range(2):
        for (b, b1) in J:
            PH = np.repeat(phi[None], grid * grid, 0); PS = np.repeat(psi[None], grid * grid, 0)
            g1, g2 = np.meshgrid(ang, ang, indexing="ij")
            PS[:, b] = g1.ravel(); PH[:, b1] = g2.ravel()
            r = I.kabsch_rmsd_batch(build_many(PH, PS), target)
            k = int(np.argmin(r))
            if r[k] < best:
                best = float(r[k]); psi[b] = PS[k, b]; phi[b1] = PH[k, b1]
    if refine and J:
        x0 = np.array([v for (b, b1) in J for v in (psi[b], phi[b1])])

        def f(x):
            ph = phi.copy(); ps = psi.copy()
            for q, (b, b1) in enumerate(J):
                ps[b] = x[2 * q]; ph[b1] = x[2 * q + 1]
            return rm(ph, ps)
        res = minimize(f, x0, method="Powell", options={"maxiter": 2000, "xtol": 1e-3, "ftol": 1e-5})
        if res.fun < best:
            best = float(res.fun)
            for q, (b, b1) in enumerate(J):
                psi[b] = res.x[2 * q]; phi[b1] = res.x[2 * q + 1]
    return best, phi, psi


def rg(W):
    W = np.asarray(W, float)
    return np.sqrt(((W - W.mean(-2, keepdims=True)) ** 2).sum(-1).mean(-1))


def split_groups(pdbs):
    f = np.array([p in I.FAIL18 for p in pdbs])
    return f, ~f


def group_means(vals, pdbs):
    vals = np.asarray(vals, float); f, o = split_groups(pdbs)
    return {"all": float(vals.mean()), "fail18": float(vals[f].mean()), "other108": float(vals[o].mean())}
