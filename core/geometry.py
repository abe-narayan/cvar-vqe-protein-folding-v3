"""Backbone geometry, superposition, PDB IO, CA-trace priors and representation floors.

Consolidates `protein_geometry.py` (628), `catrace.py` (158) and `floor.py` (73). Every
operation below has ONE implementation; where the repository carried several, the one kept
is named in its docstring together with the evidence that it is the correct one.

THE SUPERPOSITION CONVENTION, which is the thing that has actually gone wrong here.
`kabsch_superpose` builds ``R = V D U^T`` and applies it on the LEFT, as ``(R @ Pc.T).T``,
i.e. ``Pc @ R.T``. A caller that wants a matrix to POST-multiply by therefore needs
``R.T``, and getting that backwards yields a matrix that is still a valid rotation, so
nothing raises -- it just solves the downstream problem in a mis-rotated frame. That is
exactly what happened in `s10/forensics._weights_common_frame`, whose weighting solve ran
transposed and returned a "bound" worse than the equal-weight mean it is a relaxation of;
the fix moved that headline from 2.648 to 2.087. `rotation_onto` exists so the side is
stated in the signature rather than rediscovered, and `test_geometry` pins it against
`kabsch_superpose`.

THE RMSD IMPLEMENTATION KEPT. Three forms were in use: the explicit-rotation form
(`s8/audit8.rmsd_to`), the singular-value-residual form (`protein_geometry.ca_rmsd_batch`,
`s5/lib`, `s5/neighbour_ceiling`, `s7/audit`, `s8/diffuse`, `s6/refine`), and the weighted
generalisation (`s8/project.wkabsch_rmsd_batch`). Measured on all 126 cached window
universes (8.0 M window-to-native superpositions), the singular-value residual and the
explicit rotation agree to 8.9e-14 -- the earlier 4.8e-07 audit was reading a float32
cache, not a disagreement between the forms -- and the singular-value form reproduces the
cached `rr` of every one of those windows. It is kept because it never forms the rotation
(a third of the flops) and because the sign flip on the smallest singular value is the
whole reflection-rejection story in one line. `ca_rmsd_batch` here is BIT-IDENTICAL to
`protein_geometry.ca_rmsd_batch`, difference exactly 0.0.

PDB IO. `parse_pdb_ensemble` used Biopython at 26 ms/file; an exhaustive pass over
``prots/`` is 13,751 files and took 7.6 minutes of pure parsing. The column reader here is
2.4 ms/file and is validated bit-for-bit against Biopython, including the two disorder
conventions Biopython actually implements and which a naive reader gets wrong:

* a disordered ATOM resolves to the FIRST altloc of the HIGHEST occupancy (strict ``>``,
  so ties keep the earlier record) -- not the alphabetically first;
* a disordered RESIDUE -- the same sequence position deposited under two residue *names*,
  e.g. 1CBN's SER/PRO at 22 -- resolves to the LAST name encountered.

Parsed structures go through `core.cache`, keyed on the file's path, size and mtime, so a
second pass over a corpus is a load rather than a parse.
"""
from __future__ import annotations

import glob
import math
import os
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from . import cache

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BOND_N_CA = 1.458
BOND_CA_C = 1.525
BOND_C_N = 1.329
BOND_C_O = 1.231
ANGLE_N_CA_C = math.radians(111.0)
ANGLE_CA_C_N = math.radians(116.2)
ANGLE_C_N_CA = math.radians(121.7)
ANGLE_CA_C_O = math.radians(120.8)
OMEGA_TRANS = math.pi

#: Shortest N...O separation that can still be a hydrogen bond, A. Below this the donor and
#: acceptor heavy atoms are interpenetrating.
HB_MIN_ON = 2.6
#: Floor on one bond's DSSP energy. LOAD-BEARING, not redundant with HB_MIN_ON: at
#: dON = 2.6 the deepest energy this form can still reach is -5.62, in an anti-linear
#: geometry with the acceptor carbonyl C jammed 1.37 A from the donor N. -4.0 is the
#: deepest value the form assigns to a *correctly oriented* bond at the shortest
#: crystallographic N...O separation (-4.132), so the clamp costs nothing physical.
HB_E_FLOOR = -4.0

DEFAULT_PHI = math.radians(-60.0)
DEFAULT_PSI = math.radians(-45.0)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

_PDB_ACCESS_LOG: List[str] = []


def reset_pdb_log() -> None:
    _PDB_ACCESS_LOG.clear()


def get_pdb_log() -> List[str]:
    return list(_PDB_ACCESS_LOG)


# ==========================================================================
# NeRF backbone construction
# ==========================================================================
def _place_atom(a, b, c, length: float, angle: float, torsion: float):
    bcx, bcy, bcz = c[0] - b[0], c[1] - b[1], c[2] - b[2]
    nb = math.sqrt(bcx * bcx + bcy * bcy + bcz * bcz)
    if nb < 1e-9:
        nb = 1e-9
    bcx, bcy, bcz = bcx / nb, bcy / nb, bcz / nb

    abx, aby, abz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    nx = aby * bcz - abz * bcy
    ny = abz * bcx - abx * bcz
    nz = abx * bcy - aby * bcx
    nn = math.sqrt(nx * nx + ny * ny + nz * nz)
    if nn < 1e-9:
        nx, ny, nz = 0.0, 0.0, 1.0
    else:
        nx, ny, nz = nx / nn, ny / nn, nz / nn

    mx = ny * bcz - nz * bcy
    my = nz * bcx - nx * bcz
    mz = nx * bcy - ny * bcx

    d0 = -length * math.cos(angle)
    sa = length * math.sin(angle)
    d1 = sa * math.cos(torsion)
    d2 = sa * math.sin(torsion)

    return (
        c[0] + d0 * bcx + d1 * mx + d2 * nx,
        c[1] + d0 * bcy + d1 * my + d2 * ny,
        c[2] + d0 * bcz + d1 * mz + d2 * nz,
    )


def place_cb(n, ca, c):
    bx, by, bz = ca[0] - n[0], ca[1] - n[1], ca[2] - n[2]
    dx, dy, dz = c[0] - ca[0], c[1] - ca[1], c[2] - ca[2]
    ax = by * dz - bz * dy
    ay = bz * dx - bx * dz
    az = bx * dy - by * dx
    return (
        -0.58273431 * ax + 0.56802827 * bx - 0.54067466 * dx + ca[0],
        -0.58273431 * ay + 0.56802827 * by - 0.54067466 * dy + ca[1],
        -0.58273431 * az + 0.56802827 * bz - 0.54067466 * dz + ca[2],
    )


def build_backbone(phi: np.ndarray, psi: np.ndarray,
                   omega: float = OMEGA_TRANS) -> Dict[str, np.ndarray]:
    """Ideal-geometry backbone from torsions. Angstroms.

    phi[0] is never read: the chain is grown from residue 0's frame, so only psi[0..n-2]
    and phi[1..n-1] affect CA positions. psi[n-1] affects only the terminal carbonyl.
    """
    n_res = len(phi)
    if n_res < 1:
        raise ValueError("build_backbone requires at least one residue")

    N = [(0.0, 0.0, 0.0)] * n_res
    CA = [(0.0, 0.0, 0.0)] * n_res
    C = [(0.0, 0.0, 0.0)] * n_res

    N[0] = (0.0, 0.0, 0.0)
    CA[0] = (BOND_N_CA, 0.0, 0.0)
    C[0] = (CA[0][0] + BOND_CA_C * (-math.cos(ANGLE_N_CA_C)),
            CA[0][1] + BOND_CA_C * math.sin(ANGLE_N_CA_C), 0.0)

    for i in range(n_res - 1):
        N[i + 1] = _place_atom(N[i], CA[i], C[i], BOND_C_N, ANGLE_CA_C_N, psi[i])
        CA[i + 1] = _place_atom(CA[i], C[i], N[i + 1], BOND_N_CA, ANGLE_C_N_CA, omega)
        C[i + 1] = _place_atom(C[i], N[i + 1], CA[i + 1], BOND_CA_C, ANGLE_N_CA_C,
                               phi[i + 1])

    CB = [place_cb(N[i], CA[i], C[i]) for i in range(n_res)]
    O = [_place_atom(N[i], CA[i], C[i], BOND_C_O, ANGLE_CA_C_O, psi[i] + math.pi)
         for i in range(n_res)]
    return {"N": np.array(N, dtype=float), "CA": np.array(CA, dtype=float),
            "C": np.array(C, dtype=float), "CB": np.array(CB, dtype=float),
            "O": np.array(O, dtype=float)}


def _place_atom_batch(a, b, c, length: float, angle: float, torsion):
    """Vectorised `_place_atom`. Mirrors it line for line, including the degenerate-normal
    fallback to +z, so the two agree bit-for-bit rather than merely closely."""
    bc = c - b
    nb = np.linalg.norm(bc, axis=1, keepdims=True)
    np.maximum(nb, 1e-9, out=nb)
    bc = bc / nb

    ab = b - a
    nrm = np.cross(ab, bc)
    nn = np.linalg.norm(nrm, axis=1, keepdims=True)
    degenerate = (nn < 1e-9).ravel()
    nrm = nrm / np.maximum(nn, 1e-9)
    if degenerate.any():
        nrm[degenerate] = np.array([0.0, 0.0, 1.0])

    m = np.cross(nrm, bc)
    t = np.atleast_1d(np.asarray(torsion, dtype=float))
    d0 = -length * math.cos(angle)
    sa = length * math.sin(angle)
    return c + d0 * bc + (sa * np.cos(t))[:, None] * m + (sa * np.sin(t))[:, None] * nrm


def build_backbone_batch(phi: np.ndarray, psi: np.ndarray,
                         omega: float = OMEGA_TRANS) -> Dict[str, np.ndarray]:
    """`build_backbone` for a whole batch. ``phi``/``psi`` ``(B, n_res)`` -> ``(B, n, 3)``.

    The chain is sequential along the residue index but completely independent ACROSS
    structures, so the loop runs over residues while each step operates on all B at once:
    roughly ``70 * B`` Python calls become ~70 numpy operations. Numerically identical to
    the scalar builder, not approximately -- exactly.
    """
    phi = np.asarray(phi, dtype=float)
    psi = np.asarray(psi, dtype=float)
    if phi.ndim != 2 or psi.shape != phi.shape:
        raise ValueError("phi and psi must both be (B, n_res)")
    B, n_res = phi.shape
    if n_res < 1:
        raise ValueError("build_backbone_batch requires at least one residue")

    N = np.zeros((B, n_res, 3))
    CA = np.zeros((B, n_res, 3))
    C = np.zeros((B, n_res, 3))
    CA[:, 0, 0] = BOND_N_CA
    C[:, 0, 0] = CA[:, 0, 0] + BOND_CA_C * (-math.cos(ANGLE_N_CA_C))
    C[:, 0, 1] = CA[:, 0, 1] + BOND_CA_C * math.sin(ANGLE_N_CA_C)

    for i in range(n_res - 1):
        N[:, i + 1] = _place_atom_batch(N[:, i], CA[:, i], C[:, i],
                                        BOND_C_N, ANGLE_CA_C_N, psi[:, i])
        CA[:, i + 1] = _place_atom_batch(CA[:, i], C[:, i], N[:, i + 1],
                                         BOND_N_CA, ANGLE_C_N_CA, omega)
        C[:, i + 1] = _place_atom_batch(C[:, i], N[:, i + 1], CA[:, i + 1],
                                        BOND_CA_C, ANGLE_N_CA_C, phi[:, i + 1])

    b = CA - N
    d = C - CA
    CB = (-0.58273431 * np.cross(b, d) + 0.56802827 * b - 0.54067466 * d + CA)
    O = np.empty_like(CA)
    for i in range(n_res):
        O[:, i] = _place_atom_batch(N[:, i], CA[:, i], C[:, i],
                                    BOND_C_O, ANGLE_CA_C_O, psi[:, i] + math.pi)
    return {"N": N, "CA": CA, "C": C, "CB": CB, "O": O}


def amide_h_positions(N: np.ndarray, C: np.ndarray, O: np.ndarray) -> np.ndarray:
    N = np.asarray(N, dtype=float)
    C = np.asarray(C, dtype=float)
    O = np.asarray(O, dtype=float)
    n_res = len(N)
    H = np.full((n_res, 3), np.nan)
    if n_res < 2:
        return H
    d = C[:-1] - O[:-1]
    nd = np.linalg.norm(d, axis=1)
    ok = nd > 1e-6
    H[1:][ok] = N[1:][ok] + d[ok] / nd[ok][:, None]
    return H


# ==========================================================================
# Torsions
# ==========================================================================
def dihedral(p0, p1, p2, p3) -> float:
    """Scalar dihedral, radians. `dihedral_batch` is the vectorised form and is exact."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    p2, p3 = np.asarray(p2, float), np.asarray(p3, float)
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    nb1 = np.linalg.norm(b1)
    if nb1 < 1e-9:
        return 0.0
    b1n = b1 / nb1
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))


def dihedral_batch(p0, p1, p2, p3) -> np.ndarray:
    """Dihedrals of ``(..., 3)`` point stacks, radians. Degenerate ``|b1| < 1e-9`` -> 0.

    Written to mirror `dihedral` operation for operation -- same subtraction order, same
    projection, same ``atan2`` argument order -- so it reproduces it exactly rather than
    to round-off. `test_geometry.test_dihedral_batch_is_exact` asserts a difference of 0.0
    over the whole peptide database.
    """
    p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
    p2 = np.asarray(p2, float); p3 = np.asarray(p3, float)
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    nb1 = np.sqrt((b1 * b1).sum(-1))
    safe = np.where(nb1 < 1e-9, 1.0, nb1)
    b1n = b1 / safe[..., None]
    v = b0 - (b0 * b1n).sum(-1)[..., None] * b1n
    w = b2 - (b2 * b1n).sum(-1)[..., None] * b1n
    out = np.arctan2((np.cross(b1n, v) * w).sum(-1), (v * w).sum(-1))
    return np.where(nb1 < 1e-9, 0.0, out)


def extract_torsions(N, CA, C) -> Tuple[np.ndarray, np.ndarray]:
    """phi/psi from coordinates, vectorised.

    phi[0] and psi[n-1] are genuinely undefined -- they need the previous C and the next N
    -- and are filled with DEFAULT_*; encoders should exclude them rather than match
    against the placeholder.

    Was an ``n``-iteration Python loop calling `dihedral` up to twice per residue, 0.76 ms
    for a 14-mer. Two `dihedral_batch` calls, 0.05 ms, bit-identical.
    """
    N, CA, C = np.asarray(N, float), np.asarray(CA, float), np.asarray(C, float)
    n_res = len(CA)
    phi = np.zeros(n_res)
    psi = np.zeros(n_res)
    phi[0] = DEFAULT_PHI
    psi[n_res - 1] = DEFAULT_PSI
    if n_res > 1:
        phi[1:] = dihedral_batch(C[:-1], N[1:], CA[1:], C[1:])
        psi[:-1] = dihedral_batch(N[:-1], CA[:-1], C[:-1], N[1:])
    return phi, psi


def extract_torsions_batch(N, CA, C) -> Tuple[np.ndarray, np.ndarray]:
    """`extract_torsions` over a ``(B, n, 3)`` stack -> two ``(B, n)`` arrays."""
    N, CA, C = np.asarray(N, float), np.asarray(CA, float), np.asarray(C, float)
    B, n_res = N.shape[0], N.shape[1]
    phi = np.zeros((B, n_res))
    psi = np.zeros((B, n_res))
    phi[:, 0] = DEFAULT_PHI
    psi[:, n_res - 1] = DEFAULT_PSI
    if n_res > 1:
        phi[:, 1:] = dihedral_batch(C[:, :-1], N[:, 1:], CA[:, 1:], C[:, 1:])
        psi[:, :-1] = dihedral_batch(N[:, :-1], CA[:, :-1], C[:, :-1], N[:, 1:])
    return phi, psi


# ==========================================================================
# Superposition / metrics -- ONE implementation
# ==========================================================================
def _cross_covariance(Pc: np.ndarray, Qc: np.ndarray) -> np.ndarray:
    return np.einsum("bni,nj->bij", Pc, Qc)


def ca_rmsd_batch(P: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """THE canonical CA-RMSD. ``(B, n, 3)`` against one ``(n, 3)``.

    Optimal superposition without forming the rotation: with ``H = P_c^T ref_c`` and
    singular values ``s``, the residual is ``|P_c|^2 + |ref_c|^2 - 2 sum(s)``, and flipping
    the sign of the smallest singular value when ``det(V U^T) < 0`` is what forbids a
    reflection -- the step that makes a mirror image score as different rather than
    identical, which a distance matrix alone cannot do.

    Bit-identical to `protein_geometry.ca_rmsd_batch` (difference exactly 0.0) and equal to
    the explicit-rotation form of `s8/audit8.rmsd_to` to 8.9e-14, both measured over the
    8.0 M window-to-native superpositions in the 126 cached window universes.
    """
    P = np.asarray(P, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if P.ndim == 2:
        P = P[None]
    if P.ndim != 3 or reference.shape != P.shape[1:]:
        raise ValueError("P must be (B, n, 3) and reference (n, 3)")
    ref_c = reference - reference.mean(0)
    P_c = P - P.mean(1, keepdims=True)
    H = _cross_covariance(P_c, ref_c)
    U, S, Vt = np.linalg.svd(H)
    S = S.copy()
    S[:, -1] *= np.sign(np.linalg.det(np.einsum("bij,bjk->bik", U, Vt)))
    residual = (P_c ** 2).sum((1, 2)) + (ref_c ** 2).sum() - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(residual, 0.0) / P.shape[1])


def wca_rmsd_batch(P: np.ndarray, reference: np.ndarray,
                   w: np.ndarray) -> np.ndarray:
    """Per-atom-weighted CA-RMSD. Same construction, weighted centroids and covariance.

    At uniform weights it reduces to `ca_rmsd_batch` exactly, which `test_geometry` asserts;
    that equivalence is the only thing that makes a weighted "bound" comparable with the
    unweighted number it is supposed to bound.
    """
    P = np.asarray(P, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if P.ndim == 2:
        P = P[None]
    w = np.asarray(w, dtype=float)
    w = w / w.sum()
    Pc = P - (w[None, :, None] * P).sum(1, keepdims=True)
    Rc = reference - (w[:, None] * reference).sum(0, keepdims=True)
    Wp = Pc * w[None, :, None]
    H = _cross_covariance(Wp, Rc)
    U, S, Vt = np.linalg.svd(H)
    S = S.copy()
    S[:, -1] *= np.sign(np.linalg.det(np.einsum("bij,bjk->bik", U, Vt)))
    resid = (Wp * Pc).sum((1, 2)) + (w[:, None] * Rc * Rc).sum() - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(resid, 0.0))


def pairwise_ca_rmsd(P: np.ndarray, chunk: int = 128) -> np.ndarray:
    """``(B, B)`` all-pairs CA-RMSD, symmetric with a zero diagonal.

    ``B`` batched superpositions per row rather than ``B^2`` scalar ones; the chunking is
    only to bound the transient ``(chunk, B, n, 3)`` working set.
    """
    P = np.asarray(P, float)
    B = len(P)
    M = np.zeros((B, B))
    for a in range(0, B, chunk):
        for b in range(a, min(a + chunk, B)):
            M[b] = ca_rmsd_batch(P, P[b])
    return 0.5 * (M + M.T)


def kabsch_superpose(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """``mobile`` rigidly moved onto ``target``. Reflections forbidden.

    Applies its rotation on the LEFT (``(R @ Pc.T).T``). If you need a matrix to
    post-multiply by, use `rotation_onto`, not ``R``.
    """
    P = np.asarray(mobile, dtype=float)
    Q = np.asarray(target, dtype=float)
    if P.shape != Q.shape:
        raise ValueError(f"shape mismatch: {P.shape} vs {Q.shape}")
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (R @ Pc.T).T + Q.mean(axis=0)


def rotation_onto(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """The matrix ``M`` such that ``Xc @ M`` superposes centred ``X`` onto centred ``Y``.

    POST-MULTIPLY convention, stated in the signature on purpose. This is the transpose of
    the ``R`` inside `kabsch_superpose`, and using one where the other belongs is a valid
    rotation that raises nothing and silently solves everything downstream in a mis-rotated
    frame -- the defect that made `s10/forensics`' re-weighting bound worse than the
    equal-weight mean it relaxes, and moved that number from 2.648 to 2.087 once fixed.
    """
    H = np.asarray(X, float).T @ np.asarray(Y, float)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return U @ np.diag([1.0, 1.0, d]) @ Vt


def kabsch_align_batch(P: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """``(B, n, 3)`` rigidly moved onto ``(n, 3)``, in the reference's frame."""
    P = np.asarray(P, float)
    if P.ndim == 2:
        P = P[None]
    R0 = np.asarray(reference, float)
    r0 = R0.mean(0, keepdims=True)
    X, Y = P - P.mean(1, keepdims=True), R0 - r0
    U, _, Vt = np.linalg.svd(_cross_covariance(X, Y))
    d = np.sign(np.linalg.det(np.einsum("bij,bjk->bik", U, Vt)))
    D = np.zeros((len(P), 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    M = np.einsum("bij,bjk,bkl->bil", U, D, Vt)      # post-multiply, cf. rotation_onto
    return np.einsum("bni,bij->bnj", X, M) + r0


def kabsch_superpose_with_scale(mobile: np.ndarray,
                                target: np.ndarray) -> Tuple[np.ndarray, float]:
    """Superpose allowing the OPTIMAL least-squares (Umeyama) scale.

    Kept only for auditing historical numbers: now that the lattice representation returns
    Angstroms there is no legitimate reason to allow a free scale, and a similarity
    transform flatters the lattice and makes its RMSD incomparable with the torsion arm's.
    """
    P = np.asarray(mobile, dtype=float)
    Q = np.asarray(target, dtype=float)
    if P.shape != Q.shape:
        raise ValueError(f"shape mismatch: {P.shape} vs {Q.shape}")
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    var = float((Pc ** 2).sum())
    scale = float((S[:2].sum() + d * S[2]) / var) if var > 1e-12 else 1.0
    return kabsch_superpose(P * scale, Q), scale


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD of two already-superposed point sets."""
    A = np.asarray(a, dtype=float)
    B = np.asarray(b, dtype=float)
    if A.shape != B.shape:
        raise ValueError(f"shape mismatch: {A.shape} vs {B.shape}")
    return float(np.sqrt(np.mean(np.sum((A - B) ** 2, axis=1))))


def ca_rmsd(pred_ca: np.ndarray, native_ca: np.ndarray,
            allow_scale: bool = False) -> float:
    if allow_scale:
        aligned, _ = kabsch_superpose_with_scale(pred_ca, native_ca)
        return rmsd(aligned, native_ca)
    return float(ca_rmsd_batch(np.asarray(pred_ca, float)[None],
                               np.asarray(native_ca, float))[0])


def ca_rmsd_to_ensemble(pred_ca: np.ndarray,
                        ensemble_ca: Sequence[np.ndarray]) -> Dict[str, float]:
    """CA-RMSD against every deposited model of an NMR ensemble.

    For a solution structure the deposited models are all equally valid, so the min is the
    right headline and the spread across models is the floor below which a prediction
    cannot be meaningfully resolved.
    """
    vals = [ca_rmsd(pred_ca[:len(m)], np.asarray(m, float)[:len(pred_ca)])
            for m in ensemble_ca]
    arr = np.asarray(vals, dtype=float)
    return {"min": float(arr.min()), "mean": float(arr.mean()), "max": float(arr.max()),
            "model1": float(arr[0]), "best_model": int(np.argmin(arr)),
            "n_models": int(arr.size)}


def ensemble_spread(ensemble_ca: Sequence[np.ndarray]) -> float:
    """Mean pairwise CA-RMSD among deposited models -- the experimental resolution floor."""
    n = len(ensemble_ca)
    if n < 2:
        return 0.0
    stack = np.asarray([np.asarray(m, float) for m in ensemble_ca])
    M = pairwise_ca_rmsd(stack)
    iu = np.triu_indices(n, 1)
    return float(M[iu].mean())


def radius_of_gyration(coords: np.ndarray) -> float:
    c = np.asarray(coords, dtype=float)
    if c.ndim == 3:
        return np.sqrt(((c - c.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
    return float(np.sqrt(np.mean(np.sum((c - c.mean(axis=0)) ** 2, axis=1))))


def pair_distances(coords: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """``(B, npairs)`` (or ``(npairs,)``) distances for a fixed pair index set.

    The one place pairwise CA distances are computed. Indexing the pairs directly costs
    ``npairs`` rather than the ``n^2`` of a full matrix, and returns them in the flat order
    every prior, distogram and score in this repo consumes.
    """
    c = np.asarray(coords, float)
    if c.ndim == 2:
        return np.linalg.norm(c[i] - c[j], axis=-1)
    return np.linalg.norm(c[:, i, :] - c[:, j, :], axis=-1)


def distance_matrix(coords: np.ndarray) -> np.ndarray:
    """Full ``(n, n)`` or ``(B, n, n)`` Euclidean distance matrix."""
    c = np.asarray(coords, float)
    if c.ndim == 2:
        return np.linalg.norm(c[:, None, :] - c[None, :, :], axis=-1)
    return np.linalg.norm(c[:, :, None, :] - c[:, None, :, :], axis=-1)


# ==========================================================================
# Contacts and secondary structure
# ==========================================================================
def contact_map(coords: np.ndarray, threshold: float = 8.0,
                min_sep: int = 3) -> Set[Tuple[int, int]]:
    c = np.asarray(coords, dtype=float)
    n = len(c)
    if n < min_sep + 1:
        return set()
    ii, jj = np.triu_indices(n, min_sep)
    hit = pair_distances(c, ii, jj) < threshold
    return set(zip(ii[hit].tolist(), jj[hit].tolist()))


def contact_metrics(pred: Set, native: Set) -> Tuple[float, float, float]:
    if not pred and not native:
        return 1.0, 1.0, 1.0
    inter = pred & native
    p = len(inter) / len(pred) if pred else 0.0
    r = len(inter) / len(native) if native else 0.0
    f1 = 0.0 if (p + r) == 0 else 2 * p * r / (p + r)
    return p, r, f1


def dssp_energy_matrix(coords: Dict[str, np.ndarray], min_sep: int = 2,
                       cutoff: float = -0.5) -> Tuple[np.ndarray, np.ndarray]:
    """``(E, ok)``: the DSSP electrostatic energy for every ordered (donor, acceptor) pair,
    and the mask of pairs admissible as hydrogen bonds.

    THE SINGLE SOURCE OF TRUTH FOR THE DSSP FORM. It used to be written out twice -- here
    and in `energy_terms.hbond_terms` -- which is not a hypothetical drift risk: the N...O
    interpenetration fix was applied to one copy first, and until the second was found,
    structures whose donor and acceptor heavy atoms were interpenetrating were rejected by
    the energy while still being *reported* as helix or strand.
    """
    N, C, O = coords["N"], coords["C"], coords["O"]
    H = amide_h_positions(N, C, O)
    n = len(N)
    valid = np.isfinite(H).all(axis=1)

    dON = np.linalg.norm(N[:, None, :] - O[None, :, :], axis=2)
    dCH = np.linalg.norm(C[None, :, :] - H[:, None, :], axis=2)
    dOH = np.linalg.norm(H[:, None, :] - O[None, :, :], axis=2)
    dCN = np.linalg.norm(N[:, None, :] - C[None, :, :], axis=2)

    with np.errstate(divide="ignore", invalid="ignore"):
        E = 0.084 * 332.0 * (1.0 / dON + 1.0 / dCH - 1.0 / dOH - 1.0 / dCN)
    E = np.maximum(E, HB_E_FLOOR)

    idx = np.arange(n)
    sep = np.abs(idx[:, None] - idx[None, :])
    ok = valid[:, None] & (sep >= min_sep)
    ok &= (dON > HB_MIN_ON) & (dCH > 0.5) & (dOH > 0.5) & (dCN > 0.5)
    ok &= np.isfinite(E) & (E < cutoff)
    return E, ok


def dssp_hbonds(coords: Dict[str, np.ndarray], min_sep: int = 2,
                cutoff: float = -0.5) -> List[Tuple[int, int, float]]:
    E, ok = dssp_energy_matrix(coords, min_sep=min_sep, cutoff=cutoff)
    di, aj = np.where(ok)
    return [(int(i), int(j), float(E[i, j])) for i, j in zip(di, aj)]


def assign_secondary_structure(coords: Dict[str, np.ndarray]) -> str:
    """Simplified DSSP: H (helix), E (strand), C (coil).

    Helix requires *consecutive* n-turns (i->i+4 or i->i+3 at both i and i+1), which is
    what DSSP actually demands; a single i->i+4 bond over-called H.
    """
    n = len(coords["CA"])
    bonds = dssp_hbonds(coords, min_sep=2)
    ss = ["C"] * n
    nturn = np.zeros(n, dtype=bool)
    for i, j, _ in bonds:
        if j - i in (3, 4):
            nturn[i] = True
        if i - j in (3, 4):
            nturn[j] = True
    for i in range(n - 1):
        if nturn[i] and nturn[i + 1]:
            span = 5 if i + 5 <= n else n - i
            for k in range(i, min(n, i + span)):
                ss[k] = "H"
    for i, j, _ in bonds:
        if abs(i - j) >= 5 and ss[i] != "H" and ss[j] != "H":
            ss[i] = "E"
            ss[j] = "E"
    return "".join(ss)


def ss_agreement(pred: str, native: str) -> float:
    n = min(len(pred), len(native))
    if n == 0:
        return 0.0
    a = np.frombuffer(pred[:n].encode(), dtype="S1")
    b = np.frombuffer(native[:n].encode(), dtype="S1")
    return float((a == b).mean())


# ==========================================================================
# PDB IO -- one fast reader, validated against Biopython, cached on disk
# ==========================================================================
_ATOMS = ("N", "CA", "C")
_F32 = np.float32


def _read_models(path: str) -> List[Tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    """Every model of a PDB as ``(sequence, N, CA, C)``, by column-sliced ATOM parsing.

    Reproduces `Bio.PDB.PDBParser` bit-for-bit on this repository's corpora, including the
    two disorder conventions Biopython implements (see the module docstring) and the
    float32 storage of deposited coordinates -- coordinates are round-tripped through
    ``np.float32`` deliberately, because Biopython stores them that way and parsing
    ``"12.345"`` straight to float64 gives a different number.
    """
    models: List[Tuple[list, dict]] = []
    order: list = []
    cur: dict = {}
    chain_first: List[Optional[str]] = [None]

    def flush():
        if order:
            models.append((list(order), dict(cur)))
        cur.clear()
        del order[:]
        chain_first[0] = None

    with open(path, "rb") as fh:
        for raw in fh:
            rec = raw[:6]
            if rec == b"ATOM  ":
                res = raw[17:20].decode("latin-1").strip().upper()
                if res not in THREE_TO_ONE:
                    continue
                name = raw[12:16].decode("latin-1").strip()
                if name not in _ATOMS:
                    continue
                ch = raw[21:22].decode("latin-1")
                if chain_first[0] is None:
                    chain_first[0] = ch
                elif ch != chain_first[0]:
                    continue
                k = raw[22:27].decode("latin-1")          # resseq + insertion code
                d = cur.get(k)
                if d is None:
                    d = cur[k] = {"_names": []}
                    order.append(k)
                if res not in d:
                    d[res] = {}
                    d["_names"].append(res)               # DisorderedResidue: last wins
                try:
                    occ = float(raw[54:60])
                except ValueError:                        # pragma: no cover
                    occ = 0.0
                slot = d[res]
                prev = slot.get(name)
                # DisorderedAtom: strict `>` on occupancy, so ties keep the earlier record.
                if prev is not None and occ <= prev[0]:
                    continue
                slot[name] = (occ, (float(_F32(raw[30:38])), float(_F32(raw[38:46])),
                                    float(_F32(raw[46:54]))))
            elif rec == b"ENDMDL":
                flush()
    flush()

    out = []
    for order_, cur_ in models:
        seq, N, CA, C = [], [], [], []
        for k in order_:
            d = cur_[k]
            res = d["_names"][-1]
            slot = d[res]
            if not all(a in slot for a in _ATOMS):
                continue
            seq.append(THREE_TO_ONE[res])
            N.append(slot["N"][1])
            CA.append(slot["CA"][1])
            C.append(slot["C"][1])
        if seq:
            out.append(("".join(seq), np.array(N, float), np.array(CA, float),
                        np.array(C, float)))
    return out


def _read_models_biopython(path: str, chain_id: Optional[str] = None):
    """Reference reader. Only used by the equivalence test and by `parse_pdb_ensemble`
    when a specific ``chain_id`` is asked for, which the fast reader does not implement."""
    from Bio.PDB import PDBParser
    structure = PDBParser(QUIET=True).get_structure("protein", path)
    out = []
    for model in structure:
        for chain in model:
            if chain_id is not None and chain.id != chain_id:
                continue
            res = []
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                name = residue.resname.strip().upper()
                if name not in THREE_TO_ONE:
                    continue
                if not all(a in residue for a in _ATOMS):
                    continue
                res.append((THREE_TO_ONE[name], tuple(residue["N"].coord),
                            tuple(residue["CA"].coord), tuple(residue["C"].coord)))
            if not res:
                continue
            out.append(("".join(r[0] for r in res),
                        np.array([r[1] for r in res], float),
                        np.array([r[2] for r in res], float),
                        np.array([r[3] for r in res], float)))
            break
    return out


_CACHE_VERSION = 3


def parse_pdb_ensemble(path: str, chain_id: Optional[str] = None,
                       use_cache: bool = True):
    """Parse every model in a PDB -> list of ``(sequence, N, CA, C)``.

    Most short-peptide targets are solution NMR ensembles, where each deposited model is
    equally valid. Models whose sequence differs from the first are dropped rather than
    silently mixed -- the same rule Biopython's caller applied.

    Every call is recorded in the PDB access log so validation can prove the optimiser
    never touches native structures.
    """
    _PDB_ACCESS_LOG.append(os.path.abspath(path))
    if chain_id is not None:
        models = _read_models_biopython(path, chain_id=chain_id)
    else:
        st = os.stat(path)
        params = dict(path=os.path.abspath(path).replace("\\", "/").lower(),
                      size=int(st.st_size), mtime=repr(float(st.st_mtime)))
        if use_cache:
            got = cache.load("pdb", _CACHE_VERSION, **params)
            if got is not None:
                k = int(got["n_models"])
                return [(str(got["seqs"][m]), got[f"N{m}"], got[f"CA{m}"], got[f"C{m}"])
                        for m in range(k)]
        models = _read_models(path)

    kept = []
    ref_seq = None
    for seq, N, CA, C in models:
        if ref_seq is None:
            ref_seq = seq
        elif seq != ref_seq:
            break
        kept.append((seq, N, CA, C))
    if not kept:
        raise ValueError(f"No standard N/CA/C residues found in {path}")

    if chain_id is None and use_cache:
        arrays: Dict[str, np.ndarray] = {"n_models": np.array(len(kept)),
                                         "seqs": np.array([m[0] for m in kept])}
        for m, (_s, N, CA, C) in enumerate(kept):
            arrays[f"N{m}"], arrays[f"CA{m}"], arrays[f"C{m}"] = N, CA, C
        cache.store("pdb", _CACHE_VERSION, arrays, **params)
    return kept


def parse_pdb(path: str, chain_id: Optional[str] = None, model_index: int = 0):
    """One model (default the first) -> ``(sequence, N, CA, C)``."""
    models = parse_pdb_ensemble(path, chain_id=chain_id)
    return models[min(model_index, len(models) - 1)]


def _coords_from_backbone(seq, N, CA, C):
    n = len(CA)
    CB = np.array([place_cb(N[i], CA[i], C[i]) for i in range(n)])
    phi, psi = extract_torsions(N, CA, C)
    O = np.array([_place_atom(N[i], CA[i], C[i], BOND_C_O, ANGLE_CA_C_O,
                              psi[i] + math.pi) for i in range(n)])
    return seq, {"N": N, "CA": CA, "C": C, "CB": CB, "O": O}, phi, psi


def native_coords_from_pdb(path: str, chain_id: Optional[str] = None,
                           model_index: int = 0):
    """Full native coordinate dict (N, CA, C, CB, O) plus sequence and torsions."""
    return _coords_from_backbone(*parse_pdb(path, chain_id=chain_id,
                                            model_index=model_index))


def native_ensemble_from_pdb(path: str, chain_id: Optional[str] = None):
    return [_coords_from_backbone(*m)
            for m in parse_pdb_ensemble(path, chain_id=chain_id)]


def write_pdb(path: str, sequence: str, coords: Dict[str, np.ndarray],
              remark: str = "VQE predicted structure") -> None:
    lines = [f"REMARK  {remark}\n"]
    serial = 1
    for i, aa in enumerate(sequence):
        res = ONE_TO_THREE.get(aa, "GLY")
        for name, el in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
            if name not in coords:
                continue
            x, y, z = (float(v) for v in coords[name][i])
            nm = (" " + name) if len(name) < 4 else name
            lines.append(f"ATOM  {serial:>5d} {nm:<4s} {res:>3s} A{i + 1:>4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
            serial += 1
    lines.append("TER\nEND\n")
    with open(path, "w") as fh:
        fh.writelines(lines)


def scan_directory(pattern: str, verbose: bool = False):
    """Parse every PDB matching ``pattern``, yielding ``(path, models)``.

    Failures are skipped, as every caller in this repository already did. The first pass
    over a corpus pays the parse; later passes are cache loads.
    """
    for k, path in enumerate(sorted(glob.glob(pattern))):
        try:
            yield path, parse_pdb_ensemble(path)
        except Exception:
            continue
        if verbose and (k + 1) % 500 == 0:
            print(f"  parsed {k + 1}", flush=True)


# ==========================================================================
# CA-trace prior -- the term that tells a structure from its mirror image
# ==========================================================================
# The search objective is dominated by a predicted CA-CA distance matrix, and a distance
# matrix is invariant under reflection: the mirror image of any candidate has an identical
# matrix and therefore an identical score. That is a blind spot in the functional form, and
# an expensive one -- 15-24% of 1LE0's candidate pool is left-handed and on two of three
# seeds the pipeline selected one (3.32 A picked where its mirror was 2.68).
#
# The fix scores the two internal coordinates of a CA trace: the pseudo-bond-angle theta,
# which is mirror-invariant, and the pseudo-torsion tau, which CHANGES SIGN under
# reflection. Their joint distribution in deposited structures is sharply bimodal and the
# left-handed mirror of a helix lands where almost no real protein sits, so a joint
# histogram consumed as -log P rejects mirrors on physics rather than a fitted constant.
NTHETA = 24          # 0..pi
NTAU = 36            # -pi..pi
SMOOTH = 1.0         # bins of Gaussian smoothing, wrapped in tau only
FLOOR = 1e-5
CATRACE_CACHE = os.path.join(BASE, "catrace_prior.npz")


def pseudo_angles(ca: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """``(B, n-2)`` pseudo-bond-angles and ``(B, n-3)`` pseudo-torsions, radians."""
    x = np.asarray(ca, float)
    if x.ndim == 2:
        x = x[None]
    v = x[:, 1:] - x[:, :-1]
    nv = np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-9)
    u = v / nv
    theta = np.arccos(np.clip(-(u[:, :-1] * u[:, 1:]).sum(-1), -1.0, 1.0))
    b1, b2, b3 = v[:, :-2], v[:, 1:-1], v[:, 2:]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m = np.cross(n1, b2 / np.maximum(np.linalg.norm(b2, axis=-1, keepdims=True), 1e-9))
    # Negated to match IUPAC: a right-handed alpha helix must give tau = +50 deg. The
    # learned table is self-consistent either way, so this is a readability fix -- but a
    # sign convention that disagrees with the literature is a trap.
    tau = -np.arctan2((m * n2).sum(-1), (n1 * n2).sum(-1))
    return theta, tau


def _bin_angles(theta: np.ndarray, tau: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ti = np.clip((theta / np.pi * NTHETA).astype(int), 0, NTHETA - 1)
    ui = np.clip(((tau + np.pi) / (2 * np.pi) * NTAU).astype(int), 0, NTAU - 1)
    return ti, ui


def _smooth_table(h: np.ndarray) -> np.ndarray:
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


def build_catrace_table(entries) -> np.ndarray:
    """``(NTHETA, NTAU)`` negative log-probability table from held-out structures."""
    h = np.zeros((NTHETA, NTAU))
    for p in entries:
        if p.n < 4:
            continue
        th, tu = pseudo_angles(p.ca[None])
        ti, ui = _bin_angles(th[0][:-1], tu[0])          # align theta_i with tau_i
        np.add.at(h, (ti, ui), 1.0)
    h = _smooth_table(h) + FLOOR * max(h.sum(), 1.0) / h.size
    return -np.log(h / h.sum())


@lru_cache(maxsize=8)
def catrace_table(exclude_seq: str = "") -> np.ndarray:
    from . import data
    if not exclude_seq and os.path.exists(CATRACE_CACHE):
        return np.load(CATRACE_CACHE)["nlp"]
    entries = data.holdout(exclude_seq) if exclude_seq else data.load()
    nlp = build_catrace_table(entries)
    if not exclude_seq:
        np.savez(CATRACE_CACHE, nlp=nlp)
    return nlp


class CATracePrior:
    """Batched ``-log P(theta, tau)`` per structure, and a mirror diagnostic."""

    def __init__(self, exclude_seq: str = ""):
        self.nlp = catrace_table(exclude_seq)

    def score(self, ca: np.ndarray) -> np.ndarray:
        """``(B,)`` mean negative log-probability over the trace. Lower is better."""
        arr = np.asarray(ca, float)
        single = arr.ndim == 2
        if single:
            arr = arr[None]
        if arr.shape[1] < 4:
            return np.zeros(len(arr))
        th, tu = pseudo_angles(arr)
        ti, ui = _bin_angles(th[:, :-1], tu)
        out = self.nlp[ti, ui].mean(1)
        return float(out[0]) if single else out

    def handedness(self, ca: np.ndarray) -> np.ndarray:
        """``(B,)`` mean sign of the pseudo-torsion in helical windows. Positive is
        right-handed; a reflected structure returns the negation, which makes this a
        direct mirror detector rather than an indirect one."""
        arr = np.asarray(ca, float)
        if arr.ndim == 2:
            arr = arr[None]
        if arr.shape[1] < 4:
            return np.zeros(len(arr))
        th, tu = pseudo_angles(arr)
        helical = th[:, :-1] < np.radians(110.0)
        s = np.where(helical, np.sign(tu), 0.0)
        return s.sum(1) / np.maximum(helical.sum(1), 1)

    def mirror_penalty(self, ca: np.ndarray) -> np.ndarray:
        """``(B,)`` score difference between the structure and its mirror. Zero would mean
        the term cannot tell them apart, so it is reported rather than asserted."""
        arr = np.asarray(ca, float)
        if arr.ndim == 2:
            arr = arr[None]
        m = arr.copy()
        m[..., 2] *= -1.0
        return self.score(m) - self.score(arr)


# ==========================================================================
# Representation floor
# ==========================================================================
# The best CA-RMSD an encoding can express for a given native: nearest-state projection,
# batched coordinate descent, then restarts. A property of the ENCODING, measured with no
# search and no energy model -- if the floor for a target is above 2 A, no optimiser can
# put that target under 2 A and the representation is the binding constraint.
def _floor_rmsd(phi: np.ndarray, psi: np.ndarray, native_ca: np.ndarray) -> np.ndarray:
    return ca_rmsd_batch(build_backbone_batch(phi, psi)["CA"], native_ca)


# ==========================================================================
# The `numerics` backend contract (core/__init__.py)
# ==========================================================================
# `s7/audit.py` carries its own batched Kabsch, its own BLOSUM62 table and its own
# sequence encoder, deliberately retyped so that a pool built by it was not the same
# object as the pool it audits. That was the right call for an audit and is the wrong
# one for a production path, where it is three more copies to drift. These are the
# canonical implementations under the names the contract asks for. `encode` and the
# BLOSUM table live in `core.data` with the alphabet they depend on; they are reached
# here lazily so the two modules do not form an import cycle.
def kabsch_rmsd_batch(P: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Contract alias for `ca_rmsd_batch`. Bit-identical: it IS `ca_rmsd_batch`."""
    return ca_rmsd_batch(P, ref)


def pair_index(n: int, min_sep: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """The pair set every prior and distogram in this repository is indexed by."""
    return np.triu_indices(n, k=min_sep)


def pair_dists(W: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """``(B, npairs)`` CA-CA distances for a fixed pair set. See `pair_distances`."""
    return pair_distances(W, i, j)


def encode(seq: str) -> np.ndarray:
    """Sequence -> integer codes in the canonical alphabet. See `core.data.ALPHABET`."""
    from . import data
    return data.encode(seq)


def __getattr__(name):
    """Lazy re-export of the alphabet-dependent constants, avoiding an import cycle."""
    if name in ("B62", "BLOSUM62"):
        from . import data
        return data.BLOSUM62
    if name in ("AA", "ALPHABET"):
        from . import data
        return data.ALPHABET
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def project_states(rep, native_phi: np.ndarray, native_psi: np.ndarray) -> np.ndarray:
    """Per-residue nearest state by angular distance, ignoring the undefined termini.

    Vectorised over states: the angular residual is ``angle(exp(i dphi))``, computed for
    every (residue, state) cell at once instead of one Python iteration per residue.
    """
    n = rep.n_residues
    dphi = np.angle(np.exp(1j * (rep._phi - np.asarray(native_phi, float)[:, None]))) ** 2
    dpsi = np.angle(np.exp(1j * (rep._psi - np.asarray(native_psi, float)[:, None]))) ** 2
    dphi[0] = 0.0                 # phi[0] is undefined
    dpsi[n - 1] = 0.0             # psi[n-1] is undefined
    return (dphi + dpsi).argmin(1)


def descend(rep, native_ca: np.ndarray, states: np.ndarray,
            sweeps: int = 12) -> Tuple[np.ndarray, float]:
    """Coordinate descent on state assignment. Every candidate for one residue is built in
    a single `build_backbone_batch` call, so a sweep costs one batched geometry build."""
    n, k = rep.n_residues, rep.n_states
    rows = np.arange(n)
    cur = np.array(states, int)
    best = float(_floor_rmsd(rep._phi[rows, cur][None], rep._psi[rows, cur][None],
                             native_ca)[0])
    for _ in range(sweeps):
        improved = False
        for i in range(n):
            cand = np.repeat(cur[None], k, axis=0)
            cand[:, i] = np.arange(k)
            r = _floor_rmsd(rep._phi[rows[None], cand], rep._psi[rows[None], cand],
                            native_ca)
            j = int(r.argmin())
            if r[j] < best - 1e-9:
                best, cur, improved = float(r[j]), cand[j], True
        if not improved:
            break
    return cur, best


def representation_floor(rep, native_ca: np.ndarray, native_phi: np.ndarray,
                         native_psi: np.ndarray, restarts: int = 6,
                         seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    n, k = rep.n_residues, rep.n_states
    s0 = project_states(rep, native_phi, native_psi)
    rows = np.arange(n)
    proj = float(_floor_rmsd(rep._phi[rows, s0][None], rep._psi[rows, s0][None],
                             native_ca)[0])
    best_s, best = descend(rep, native_ca, s0)
    dproj = best
    for _ in range(restarts):
        start = np.where(rng.random(n) < 0.35, rng.integers(0, k, n), best_s)
        s, v = descend(rep, native_ca, start)
        if v < best:
            best, best_s = v, s
    return {"projection": proj, "descent_from_projection": dproj,
            "floor": best, "states": best_s}
