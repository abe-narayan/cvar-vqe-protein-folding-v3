"""`core.project` -- STAGE 3b, the projection onto the manifold of ideal-geometry chains.

What this stage is
The coordinate average of the filtered pool is not a peptide: its mean CA-CA bond is
~2.96 A against a real 3.80.  This stage finds the torsion vector whose ideal-geometry
backbone is nearest that point cloud,

    minimise   CA-RMSD( build(phi, psi), C )  +  lam * pen(phi, psi)

over ``(phi, psi)``, by L-BFGS-B from several starts.  It is `s8.project` consolidated --
same objective, same penalty, same starts, same ladder, same L-BFGS-B call -- with the
inner loop rewritten.  `tests/test_project.py` measures the equivalence on real pipeline
inputs rather than asserting it.

Why it needed rewriting -- the measured profile
Measured on the 126-target instrument at production settings (K=500, m=75, `ramah`,
lam=0.3, maxiter=300, multi_start=True), the projection was 43.1% of the whole pipeline at
8.2 s per target, and `cProfile` on the reference put 93.4% of that inside ONE function:

    build_backbone_batch          101.2 s of 108.4 s cumulative  (93.4%)
      _place_atom_batch            786,390 calls
        numpy.cross              1,589,366 calls, 65.5 s cumulative
    kabsch_rmsd_batch               4.2 s  ( 3.8%)
    _bilinear (the penalty)         0.6 s  ( 0.6%)

The optimiser was NOT the cost and neither was the penalty.  The cost was that the
reference forward map is a Python loop over residues, three `_place_atom_batch` calls per
residue plus one for the carbonyl oxygen, each of which calls `numpy.cross` twice, and
`numpy.cross` spends most of its time in `moveaxis` / `normalize_axis_tuple` deciding what
it was asked to do.  ~47 `_place_atom_batch` calls x ~60 numpy-level calls each = ~2,800
numpy dispatches to build one 13-residue chain, and there are ~2,100 builds per target.

What was done about it
Three things make the forward map cheaper.  Only the first two are in the SHIPPED path,
because only the first two cannot change a value -- see the next section, which is the
part of this docstring that was written after an audit proved the third one does.

1.  **The carbonyl oxygen and CB are not built.**  The objective reads `CA` and nothing
    else; the reference builds `O` in a second n-step loop and `CB` from a cross product,
    26.6% of the placement calls, and throws both away.  Removing work that no value
    depends on cannot move a number.

2.  **`numpy.cross` is written out** (`_cross3`) as the permutation
    ``a[(1,2,0)] b[(2,0,1)] - a[(2,0,1)] b[(1,2,0)]`` -- the same two products and the same
    subtraction, 7 numpy dispatches instead of ~15, asserted identical at exactly 0.0.  The
    two norms lose `numpy.linalg.norm`'s dispatch the same way, as the `sqrt` of the
    `add.reduce` it performs internally.  Together with (1) this is `build_ca_exact`, which
    is bit-identical to `build_backbone_batch` at every batch width, and it is ~2.2x.

3.  **A matrix scan forward map, and an analytic gradient.**  Every bond length and bond
    angle in this builder is a constant; only the torsion varies.  So placing the next atom
    is a FIXED 4x4 rigid transform parameterised by one torsion, the chain is the cumulative
    product of those transforms, and a cumulative product is an associative scan:
    ``ceil(log2(3(n-1)))`` batched 4x4 matmuls instead of ``3(n-1)`` sequential placements.
    The same transforms give the gradient for free: a torsion rotates its downstream
    sub-chain rigidly about a bond axis that is a column of the cumulative frame, so
    ``d(CA_k)/dt = u_s x (CA_k - a_s)``, and reordering through ``g.(u x v) = u.(v x g)``
    with two suffix sums makes ALL 2n derivatives O(n) instead of 2n forward builds.
    Together: 13x, agreeing with central differences to 1.6e-9 and with the reference
    builder to 3.1e-13 A at n=40.  AND IT IS NOT THE DEFAULT.  Read on.

Three modes, and why the fast one is not the default
This is the correction that matters, and it was found by an independent audit rather than
by this module, which had claimed it "changes only how fast one L-BFGS-B iteration is
computed".  That is true of the science knobs and FALSE OF THE OUTPUT.

    ``exact``     bit-exact builder + the canonical Kabsch + the reference's own one-sided
                  finite difference at the reference's eps.  Reproduces `s8.project`
                  BIT-FOR-BIT, 0.0 on 126/126 targets.  **THE DEFAULT.**  2.2x.
    ``analytic``  scan builder + the true gradient.  13.3x, and it emits DIFFERENT
                  STRUCTURES: 126/126 targets move, median 0.042 A, worst 1.90 A.
    ``fd``        scan builder + the reference's finite difference.  Diagnostic only.

The reason is the degeneracy this module's own docstring already warned about, biting one
level lower than expected.  The scan builder agrees with the reference to ~1e-13 A, which
is nothing -- but L-BFGS-B is started from a fully extended or fully helical chain, a long
way from any minimum, and its early steps are large.  A 1e-13 A difference in the forward
map is enough to route some trajectories into the OTHER ideal-geometry torsion branch, and
the two branches sit at near-equal objective distance from the average.  Measured on the
126-target instrument with the reference's own gradient formula, so that the gradient is
not the variable: the scan builder alone moves the emitted structure on 126/126 targets,
median 0.031 A, worst 1.63 A, while the objective it reaches is lower on 66 and higher on
60 -- it is not converging worse, it is landing somewhere else.

So "faster and equivalent" was not available, and the choice is stated rather than made
quietly.  `exact` ships.  `analytic` is 4.6x faster again and is a legitimate arm to
evaluate, but it is a DIFFERENT PIPELINE and is keyed as one.

Which mode is in the cache key
`GRAD` is a module default and is NOT read from the environment.  An import-time
environment global cannot reach `Config.key()`, and the audit demonstrated the
consequence: two runs of `smoke8` under different modes produced the same `cfg_key`, and
the second "completed" all 8 targets in 0.215 s by serving the first's arrays.  The mode
travels as `core.pipeline.Config.project_grad`, is hashed like every other parameter that
changes a number, and is passed per call as `lam_path(..., grad=...)` so it cannot leak
between targets inside a worker.

The hard constraint -- the branch structure is load-bearing
The projection is DEGENERATE.  A CA trace admits two ideal-geometry torsion solutions at
near-equal objective distance, one Ramachandran-plausible and one not, and a warm-started
optimiser cannot cross between them (S9-2).  That is why `lam_path(multi=True)` puts the
four generic starts back at every rung.  **The number of starts, `maxiter`, the penalty,
`lam` and the convergence tolerance are science, not tuning knobs, and nothing here
touches them.**

And the degeneracy is not only a warning about multi-start.  `python -m core.project
stability` asks `s8.project` itself to project onto ``C`` and onto ``R C + t``, which is
the same problem -- a CA-RMSD is exactly invariant under a rigid motion and the penalty
never sees the cloud -- with every floating-point operation in the Kabsch differently
rounded.  THE REFERENCE DISAGREES WITH ITSELF, by up to 1.6 A on the same targets where
the scan builder disagrees with it.  So the structures this stage returns on those targets
are not determined by the objective; they are determined by the arithmetic.  That is a
property of the stage, it is measured rather than asserted, and it is the reason the
shipped mode is bit-exact instead of merely accurate.

Related: the distance objective is exactly mirror-blind, so a distance-based multi-start
selects enantiomers.  This objective is a COORDINATE distance and is chirality-sensitive,
which is what makes lowest-objective selection safe here.  `test_project.py` asserts every
emitted structure is L-handed anyway, because that safety is a property of the objective
and would be silently lost if the objective were ever reformulated.

    python -m core.project selfcheck    # gradient vs finite differences, builder vs ref
    python -m core.project harvest      # persist the real stage-3b input for 126 targets
    python -m core.project equiv        # the four-arm per-target comparison
    python -m core.project exactness    # is the shipped mode bit-identical on all 126?
    python -m core.project stability    # the reference against ITSELF under a null change
    python -m core.project degeneracy   # the objective gap between the best two starts
    python -m core.project iters        # is maxiter=300 binding?
    python -m core.project bench N      # time the stage on N real harvested inputs
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

import core

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(HERE)

# ============================================================ geometry constants
#: Identical to `core.geometry` / `protein_geometry`.  Bound here rather than imported so
#: the scan builder cannot silently disagree with the module it must reproduce -- the
#: self-check asserts they are the same objects' values.
BOND_N_CA = 1.458
BOND_CA_C = 1.525
BOND_C_N = 1.329
ANGLE_N_CA_C = math.radians(111.0)
ANGLE_CA_C_N = math.radians(116.2)
ANGLE_C_N_CA = math.radians(121.7)
OMEGA_TRANS = math.pi

#: The three placements per residue, in build order: psi -> N(i+1), omega -> CA(i+1),
#: phi(i+1) -> C(i+1).  `(bond angle, bond length)` for each.
_STEP_ANGLE = (ANGLE_CA_C_N, ANGLE_C_N_CA, ANGLE_N_CA_C)
_STEP_BOND = (BOND_C_N, BOND_N_CA, BOND_CA_C)

#: The gradient modes, and what each one costs in fidelity.  See `_make_fg`.
GRAD_MODES = ("exact", "analytic", "fd")

#: The DEFAULT is `exact`, which reproduces `s8.project` bit-for-bit.  It is NOT read from
#: the environment: the mode changes the emitted structure, so it has to travel in
#: `Config.key()` like every other parameter that changes a number, and an import-time
#: environment global cannot.  `core.pipeline.Config.project_grad` carries it, and
#: `SELECTABLE_GRADIENT` is how the pipeline knows this backend can be told.
GRAD = "exact"
SELECTABLE_GRADIENT = True


def set_gradient(mode: str) -> str:
    """Set the module default.  Prefer passing `grad=` to `lam_path`, which is per-call
    and cannot leak between targets."""
    global GRAD
    if mode not in GRAD_MODES:
        raise KeyError(f"unknown gradient mode {mode!r}; known: {GRAD_MODES}")
    GRAD = mode
    return GRAD

#: The reference's finite-difference step, in radians.  PINNED -- it is part of the
#: reference optimisation, not a tolerance.
FD_EPS = 1e-5

# ============================================================ the forward map
#: The two cyclic permutations a cross product is: ``(a x b)_i = a_{i+1} b_{i+2} -
#: a_{i+2} b_{i+1}``.  Bound once because building them per call costs more than the
#: multiply does.
_CYC1 = np.array([1, 2, 0])
_CYC2 = np.array([2, 0, 1])


def _cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """`numpy.cross` for ``(..., 3)`` only, written out.

    NOT a micro-optimisation for its own sake: `numpy.cross` was 60% of the reference
    builder's runtime and 31% of the optimised gradient's, and almost none of that is
    arithmetic -- it is `moveaxis` / `normalize_axis_tuple` working out which axis it was
    handed.  The products and the subtraction are in the same order `numpy.cross` uses, so
    the result is bit-identical, which `tests/test_project.py` asserts at exactly 0.0.
    """
    return a[..., _CYC1] * b[..., _CYC2] - a[..., _CYC2] * b[..., _CYC1]


def _det3(M: np.ndarray) -> float:
    """``det`` of one 3x3 by cofactor expansion.  Used only for its SIGN, on a product of
    orthogonal matrices whose determinant is exactly +-1, so it is nowhere near the
    cancellation regime where a LAPACK LU would differ."""
    return float(M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
                 - M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0])
                 + M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0]))


def _initial_frame() -> np.ndarray:
    """``(4, 4)`` frame at C[0]: columns ``[bc, m, nrm]``, translation ``C[0]``.

    The reference seeds the chain with N[0] = origin, CA[0] = (BOND_N_CA, 0, 0) and C[0]
    placed in the xy-plane, then makes its first call `_place_atom(N[0], CA[0], C[0], ...)`.
    That call's frame is ``bc = unit(C0 - CA0)``, ``nrm = unit((CA0 - N0) x bc) = +z``,
    ``m = nrm x bc``, anchored at C[0] -- which is this matrix.
    """
    ca0 = np.array([BOND_N_CA, 0.0, 0.0])
    c0 = ca0 + np.array([BOND_CA_C * (-math.cos(ANGLE_N_CA_C)),
                         BOND_CA_C * math.sin(ANGLE_N_CA_C), 0.0])
    bc = c0 - ca0
    bc = bc / np.linalg.norm(bc)
    nrm = np.array([0.0, 0.0, 1.0])
    m = np.cross(nrm, bc)
    G = np.zeros((4, 4))
    G[:3, 0] = bc
    G[:3, 1] = m
    G[:3, 2] = nrm
    G[:3, 3] = c0
    G[3, 3] = 1.0
    return G

_G0 = _initial_frame()

_GEOM_CACHE: dict = {}


def _step_geometry(M: int):
    """``(cos th, sin th, bond)`` per scan step, cached by chain length."""
    got = _GEOM_CACHE.get(M)
    if got is None:
        cth = np.empty(M)
        sth = np.empty(M)
        ln = np.empty(M)
        for k in range(3):
            cth[k::3] = math.cos(_STEP_ANGLE[k])
            sth[k::3] = math.sin(_STEP_ANGLE[k])
            ln[k::3] = _STEP_BOND[k]
        got = (cth, sth, ln)
        _GEOM_CACHE[M] = got
    return got


def _local_transforms(t: np.ndarray) -> np.ndarray:
    """``(B, M, 4, 4)`` local rigid transforms for scan torsions ``t`` ``(B, M)``.

    DERIVATION.  `_place_atom(a, b, c, l, th, t)` writes the new atom as

        d = c + F . (-l cos th,  l sin th cos t,  l sin th sin t)

    with the orthonormal frame ``F = [bc, m, nrm]``, ``bc = unit(c - b)``,
    ``nrm = unit((b - a) x bc)``, ``m = nrm x bc``.  The NEXT placement's frame is built
    from ``(b, c, d) = (b', c', d')``, and because ``|d - c| = l`` and ``a' - b'`` is
    parallel to ``bc``, every one of its columns is a FIXED linear combination of this
    frame's columns:

        bc' = F u,   u = (-cos th,  sin th cos t,  sin th sin t)
        nrm' = F q,  q = ( 0,      -sin t,         cos t       )
        m'  = F w,   w = (-sin th, -cos th cos t, -cos th sin t)
        c'  = c + F (l u)

    so ``G' = G T`` with ``T = [[u | w | q | l u], [0 0 0 1]]``, and the whole chain is the
    cumulative product ``G_s = G_0 T_1 ... T_s`` whose translation column is atom ``s``.
    ``T`` is a rotation about the frame's own ``bc`` axis composed with a fixed screw, and
    it depends on the conformation ONLY through the single torsion ``t``.
    """
    B, M = t.shape
    cth, sth, ln = _step_geometry(M)
    ct = np.cos(t)
    st = np.sin(t)
    T = np.zeros((B, M, 4, 4))
    sc = sth * ct
    ss = sth * st
    T[:, :, 0, 0] = -cth
    T[:, :, 1, 0] = sc
    T[:, :, 2, 0] = ss
    T[:, :, 0, 1] = -sth
    T[:, :, 1, 1] = -cth * ct
    T[:, :, 2, 1] = -cth * st
    T[:, :, 1, 2] = -st
    T[:, :, 2, 2] = ct
    T[:, :, 0, 3] = ln * (-cth)
    T[:, :, 1, 3] = ln * sc
    T[:, :, 2, 3] = ln * ss
    T[:, :, 3, 3] = 1.0
    return T


def frames(phi: np.ndarray, psi: np.ndarray, omega: float = OMEGA_TRANS) -> np.ndarray:
    """``(B, M+1, 4, 4)`` cumulative frames, ``M = 3(n-1)``.

    ``G[:, s, :3, 3]`` is atom ``s`` of the flat chain ``C0, N1, CA1, C1, N2, ...`` and
    ``G[:, s, :3, 0]`` is the axis of the torsion applied at step ``s+1``, anchored at that
    same atom.  Both are what the analytic gradient needs, so the gradient costs no extra
    forward work.

    The scan is Hillis-Steele: ``log2(M)`` rounds of ``P[d:] = P[:M-d] @ P[d:]``, which
    keeps the left-to-right order matrix products require.
    """
    phi = np.asarray(phi, float)
    psi = np.asarray(psi, float)
    if phi.ndim != 2 or psi.shape != phi.shape:
        raise ValueError("phi and psi must both be (B, n_res)")
    B, n = phi.shape
    if n < 1:
        raise ValueError("frames requires at least one residue")
    M = 3 * (n - 1)
    if M == 0:
        return np.broadcast_to(_G0, (B, 1, 4, 4)).copy()
    t = np.empty((B, M))
    t[:, 0::3] = psi[:, :n - 1]
    t[:, 1::3] = omega
    t[:, 2::3] = phi[:, 1:]
    P = _local_transforms(t)
    buf = np.empty_like(P)
    d = 1
    while d < M:
        np.matmul(P[:, :M - d], P[:, d:], out=buf[:, d:])
        buf[:, :d] = P[:, :d]
        P, buf = buf, P
        d *= 2
    G = np.empty((B, M + 1, 4, 4))
    G[:, 0] = _G0
    np.matmul(_G0, P, out=G[:, 1:])
    return G


def _place_exact(a, b, c, length: float, angle: float, torsion):
    """`core.geometry._place_atom_batch`, BIT-FOR-BIT, with the dispatch taken out.

    Same operations in the same order on the same values -- the only changes are
    `numpy.cross` written out (`_cross3`, asserted identical at 0.0), the two norms written
    as the `sqrt` of an `add.reduce` that `numpy.linalg.norm` performs anyway, and the
    degeneracy branch tested with one `min` instead of a comparison, a `ravel` and an `any`.
    Nothing here is an approximation of the reference; it IS the reference, dispatched less.
    """
    bc = c - b
    nb = np.sqrt((bc * bc).sum(1))[:, None]
    np.maximum(nb, 1e-9, out=nb)
    bc = bc / nb

    nrm = _cross3(b - a, bc)
    nn = np.sqrt((nrm * nrm).sum(1))[:, None]
    if nn.min() < 1e-9:                                 # never fires on a real chain
        deg = (nn < 1e-9).ravel()
        nrm = nrm / np.maximum(nn, 1e-9)
        nrm[deg] = np.array([0.0, 0.0, 1.0])
    else:
        nrm = nrm / nn

    m = _cross3(nrm, bc)
    d0 = -length * math.cos(angle)
    sa = length * math.sin(angle)
    return (c + d0 * bc + (sa * np.cos(torsion))[:, None] * m
            + (sa * np.sin(torsion))[:, None] * nrm)


def build_ca_exact(phi: np.ndarray, psi: np.ndarray, omega: float = OMEGA_TRANS):
    """``(B, n, 3)`` CA trace, BIT-IDENTICAL to `build_backbone_batch(...)["CA"]`.

    The point of this function is the word identical.  The scan builder is 20x faster and
    agrees only to ~1e-13 A, and that is not good enough here: the projection is degenerate,
    L-BFGS-B started from a fully extended or fully helical chain takes large early steps,
    and a 1e-13 A difference in the forward map routes some trajectories into the other
    torsion branch.  Measured on the 126-target instrument, the scan builder under the
    reference's own finite-difference gradient moves the emitted structure on 126/126
    targets, median 0.031 A and worst 1.63 A.  So the fast builder cannot be the default,
    and this one -- which reproduces the reference exactly -- is.

    The speed here comes only from what does not change a value: `numpy.cross` written out,
    the carbonyl oxygen and CB not built at all (the objective reads CA and the reference
    threw both away), and the norms not routed through `numpy.linalg.norm`'s dispatch.

    Where this stops, and why it is not physics.  `_place_exact` is now about 60 numpy
    dispatches averaging ~0.47 us each, which is the dispatch cost itself: 85% of the
    shipped objective is this function, and going faster means issuing fewer numpy calls,
    which means changing the sequence of operations, which is the one thing that costs
    bit-identity.  IN NUMPY that is the floor.  It is NOT the floor in general -- a
    compiled loop performing the SAME operations in the SAME order would be bit-identical,
    because IEEE754 is deterministic for a fixed operation sequence, and it would be far
    faster.  That route is blocked on tooling rather than closed on principle: numba and
    cython are not installed, no C compiler is on PATH, and installing numba is actively
    the wrong move -- it carries a hard numpy upper bound against this environment's numpy
    2.5.1, so a resolver would likely DOWNGRADE numpy and change summation and BLAS
    behaviour underneath every pinned golden in the repository.  Recorded here so that if a
    compiler ever appears, the next person knows this is open.
    """
    phi = np.asarray(phi, dtype=float)
    psi = np.asarray(psi, dtype=float)
    if phi.ndim != 2 or psi.shape != phi.shape:
        raise ValueError("phi and psi must both be (B, n_res)")
    B, n = phi.shape
    if n < 1:
        raise ValueError("build_ca_exact requires at least one residue")
    N = np.zeros((B, n, 3))
    CA = np.zeros((B, n, 3))
    C = np.zeros((B, n, 3))
    CA[:, 0, 0] = BOND_N_CA
    C[:, 0, 0] = CA[:, 0, 0] + BOND_CA_C * (-math.cos(ANGLE_N_CA_C))
    C[:, 0, 1] = CA[:, 0, 1] + BOND_CA_C * math.sin(ANGLE_N_CA_C)
    om = np.full(B, omega)
    for i in range(n - 1):
        N[:, i + 1] = _place_exact(N[:, i], CA[:, i], C[:, i],
                                   BOND_C_N, ANGLE_CA_C_N, psi[:, i])
        CA[:, i + 1] = _place_exact(CA[:, i], C[:, i], N[:, i + 1],
                                    BOND_N_CA, ANGLE_C_N_CA, om)
        C[:, i + 1] = _place_exact(C[:, i], N[:, i + 1], CA[:, i + 1],
                                   BOND_CA_C, ANGLE_N_CA_C, phi[:, i + 1])
    return CA


def build_ca(phi: np.ndarray, psi: np.ndarray, omega: float = OMEGA_TRANS):
    """``(CA, G)`` -- the CA trace ``(B, n, 3)`` and the frames behind it.

    Reproduces `core.geometry.build_backbone_batch(...)["CA"]` to ~1e-13 A; see the module
    docstring on why nothing emitted comes from here.
    """
    phi = np.asarray(phi, float)
    B, n = phi.shape
    G = frames(phi, psi, omega)
    CA = np.empty((B, n, 3))
    CA[:, 0, 0] = BOND_N_CA
    CA[:, 0, 1:] = 0.0
    if n > 1:
        CA[:, 1:] = G[:, 2::3, :3, 3]
    return CA, G


# ============================================================ the data term
def _rmsd_batch(CA: np.ndarray, ref_c: np.ndarray, ref_sq: float) -> np.ndarray:
    """CA-RMSD of every ``(n, 3)`` in `CA` against a PRE-CENTRED reference.

    The construction `s7.audit.kabsch_rmsd_batch` and `core.geometry.ca_rmsd_batch` both
    use: residual = |P_c|^2 + |ref_c|^2 - 2 sum(s), with the smallest singular value
    sign-flipped when det(V U^T) < 0, which is the step that forbids a reflection and is
    the whole reason this objective can tell an L chain from a D one.
    """
    Pc = CA - CA.mean(1, keepdims=True)
    H = np.einsum("bni,nj->bij", Pc, ref_c)
    U, S, Vt = np.linalg.svd(H)
    S = S.copy()
    S[:, -1] *= np.sign(np.linalg.det(np.einsum("bji,bkj->bik", Vt, U)))
    resid = (Pc ** 2).sum((1, 2)) + ref_sq - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(resid, 0.0) / CA.shape[1])


def _rmsd_and_dca(CA1: np.ndarray, ref_c: np.ndarray, ref_sq: float):
    """``(rmsd, d rmsd / d CA)`` for ONE structure ``(n, 3)``.

    The superposition is optimal, so by the envelope theorem the derivative through the
    rotation and the translation both vanish and

        d rmsd / d p_j = (p_j,c - R^T r_j,c) / (n rmsd)

    with ``R = V diag(1, 1, det) U^T`` the same rotation whose singular values the residual
    was read off.  The centroid term drops out exactly because both sets are centred:
    ``sum_i (R p_i,c - r_i,c) = 0``.
    """
    n = CA1.shape[0]
    Pc = CA1 - CA1.mean(0)
    H = Pc.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    dsign = 1.0 if _det3(V @ U.T) >= 0.0 else -1.0
    resid = (Pc ** 2).sum() + ref_sq - 2.0 * (S[0] + S[1] + dsign * S[2])
    rmsd = math.sqrt(max(resid, 0.0) / n)
    if rmsd <= 0.0:
        return 0.0, np.zeros_like(Pc)
    Vd = V.copy()
    Vd[:, 2] *= dsign
    R = Vd @ U.T                                      # p -> R p aligns onto ref
    return rmsd, (Pc - ref_c @ R) / (n * rmsd)


def _torsion_grad(G: np.ndarray, CA1: np.ndarray, g: np.ndarray) -> np.ndarray:
    """``(2n,)`` d(objective)/d(phi, psi) from ``g = d(objective)/d(CA)``.

    A torsion applied at scan step ``s+1`` rotates every atom downstream of it rigidly
    about the axis ``u_s`` through the anchor ``a_s`` -- both columns of the cumulative
    frame ``G[s]``.  So ``d CA_k / d t = u_s x (CA_k - a_s)`` for the affected suffix, and

        d f / d t = sum_k g_k . (u_s x (CA_k - a_s))
                  = u_s . ( sum_k CA_k x g_k  -  a_s x sum_k g_k )

    by ``g . (u x v) = u . (v x g)``.  The two sums are suffix sums over the CA index, so
    all 2n derivatives cost one reverse cumulative sum and two dot products -- O(n), not
    O(n) forward builds.

    ``phi[0]`` and ``psi[n-1]`` are never read by the builder, so their entries are exactly
    zero, which is also exactly what the reference's one-sided difference returns for them.
    """
    n = CA1.shape[0]
    out = np.zeros(2 * n)
    if n < 2:
        return out
    Ssuf = np.zeros((n + 1, 3))
    Qsuf = np.zeros((n + 1, 3))
    np.cumsum(g[::-1], axis=0, out=Ssuf[:n][::-1])
    np.cumsum(_cross3(CA1, g)[::-1], axis=0, out=Qsuf[:n][::-1])

    #: psi[i] (i = 0..n-2) is applied at step 3i+1, whose frame is G[3i]; it moves CA[i+1:]
    Gp = G[0:3 * (n - 1):3]
    Vp = Qsuf[1:n] - _cross3(Gp[:, :3, 3], Ssuf[1:n])
    out[n:2 * n - 1] = np.einsum("ij,ij->i", Gp[:, :3, 0], Vp)

    #: phi[j] (j = 1..n-1) is applied at step 3j, whose frame is G[3j-1]; it moves CA[j+1:]
    Gf = G[2:3 * (n - 1) + 1:3]
    Vf = Qsuf[2:n + 1] - _cross3(Gf[:, :3, 3], Ssuf[2:n + 1])
    out[1:n] = np.einsum("ij,ij->i", Gf[:, :3, 0], Vf)
    return out

# ============================================================ the torsion priors
#: Ramachandran grid, smoothing width in bins, and the Dirichlet pseudo-count.  Copied by
#: VALUE from `s8.project` -- these are part of the trained prior, not parameters.
RB = 36
SIGMA_BINS = 1.5
PSEUDO = 0.5
HINGE_PCT = 5.0
PHIP_TAU = 0.15
CLASSES = ("gen", "gly", "pro", "prepro")
AA20 = "ARNDCQEGHILKMFPSTWYV"
PENALTIES = ("rama", "rama20", "vm", "phip", "ramah")

#: The fold-disciplined count tables, accumulated once by `s8.project` from
#: `peptide_db.load()` restricted to the training folds plus `distogram._fold_fragments`.
#: They are DATA, like the distogram checkpoints, and are read from where that module
#: wrote them; if the file is absent it is rebuilt through the same code path.
PRIOR_JSON = os.path.join(_ROOT, "s8", "project_prior.json")


def res_classes(seq) -> np.ndarray:
    """``(n,)`` GLY / PRO / PRE-PRO / GENERAL class index.  Priority order as declared:
    a glycine followed by proline scores as glycine."""
    n = len(seq)
    out = np.zeros(n, int)
    for i, c in enumerate(seq):
        if c == "G":
            out[i] = 1
        elif c == "P":
            out[i] = 2
        elif i + 1 < n and seq[i + 1] == "P":
            out[i] = 3
    return out


def _raw_tables() -> dict:
    if not os.path.exists(PRIOR_JSON):                              # pragma: no cover
        from s8 import project as _s8
        return _s8.build_tables(verbose=False)
    with open(PRIOR_JSON) as fh:
        return json.load(fh)


def _wrap_smooth(C: np.ndarray, sigma_bins: float = SIGMA_BINS) -> np.ndarray:
    """Wrapped-Gaussian smoothing of a ``(..., RB, RB)`` count table, by FFT."""
    k = np.arange(RB)
    d = np.minimum(k, RB - k)
    g = np.exp(-0.5 * (d / sigma_bins) ** 2)
    g = g / g.sum()
    G = np.fft.rfft(g)
    A = np.fft.irfft(np.fft.rfft(C, axis=-1) * G, n=RB, axis=-1)
    A = np.fft.irfft(np.fft.rfft(A, axis=-2) * G[:, None], n=RB, axis=-2)
    return np.maximum(A, 0.0)

_TAB: dict = {}
_HINGE: dict = {}


def logp_tables(kind: str = "rama") -> np.ndarray:
    """``(5, C, RB, RB)`` log-densities: smoothed, pseudo-counted, normalised per class."""
    if kind not in _TAB:
        raw = _raw_tables()
        C = np.asarray(raw["cls" if kind == "rama" else "aa"], float)
        S = _wrap_smooth(C) + PSEUDO
        _TAB[kind] = np.log(S / S.sum((2, 3), keepdims=True))
    return _TAB[kind]


def hinge_thresholds(pct: float = HINGE_PCT) -> np.ndarray:
    """``(5, C)`` per-fold, per-class log-density below which a residue is implausible:
    the `pct`-th percentile over REAL RESIDUES of that class, not over grid area."""
    key = round(pct, 6)
    if key not in _HINGE:
        cnt = np.asarray(_raw_tables()["cls"], float)
        L = logp_tables("rama")
        out = np.zeros(L.shape[:2])
        for f in range(L.shape[0]):
            for k in range(L.shape[1]):
                v = L[f, k].ravel()
                w = cnt[f, k].ravel()
                o = np.argsort(v)
                cw = np.cumsum(w[o])
                tot = cw[-1]
                if tot <= 0:
                    out[f, k] = v.min()
                    continue
                out[f, k] = v[o[int(np.searchsorted(cw, tot * pct / 100.0))]]
        _HINGE[key] = out
    return _HINGE[key]


def _bilinear(L: np.ndarray, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
    """Bilinear interpolation of a per-residue ``(n, RB, RB)`` LOG-density.  Wrapped
    exactly on both axes; interpolating the log keeps the penalty finite everywhere."""
    u = (np.asarray(phi, float) + math.pi) / (2 * math.pi) * RB - 0.5
    v = (np.asarray(psi, float) + math.pi) / (2 * math.pi) * RB - 0.5
    i0 = np.floor(u).astype(int)
    j0 = np.floor(v).astype(int)
    fu, fv = u - i0, v - j0
    i0 %= RB
    j0 %= RB
    i1, j1 = (i0 + 1) % RB, (j0 + 1) % RB
    r = np.arange(L.shape[0])[None, :]
    return ((1 - fu) * (1 - fv) * L[r, i0, j0] + fu * (1 - fv) * L[r, i1, j0]
            + (1 - fu) * fv * L[r, i0, j1] + fu * fv * L[r, i1, j1])

_ROWS: dict = {}


def _rows(k: int) -> np.ndarray:
    """``arange(k)[None, :]``, cached -- the per-residue table index, rebuilt on every
    penalty evaluation otherwise."""
    got = _ROWS.get(k)
    if got is None:
        got = _ROWS[k] = np.arange(k)[None, :]
    return got


def _bilinear_d(L: np.ndarray, phi: np.ndarray, psi: np.ndarray):
    """``(lp, dlp/dphi, dlp/dpsi)``.  The interpolant is piecewise bilinear, so its
    derivative is exact inside a cell and one-sided on a cell boundary -- which is where
    it differs, by a bounded amount, from a finite difference that straddles the seam."""
    s = RB / (2 * math.pi)
    u = (np.asarray(phi, float) + math.pi) * s - 0.5
    v = (np.asarray(psi, float) + math.pi) * s - 0.5
    i0 = np.floor(u).astype(int)
    j0 = np.floor(v).astype(int)
    fu, fv = u - i0, v - j0
    i0 %= RB
    j0 %= RB
    i1, j1 = (i0 + 1) % RB, (j0 + 1) % RB
    r = _rows(L.shape[0])
    a = L[r, i0, j0]
    b = L[r, i1, j0]
    c = L[r, i0, j1]
    d = L[r, i1, j1]
    lp = (1 - fu) * (1 - fv) * a + fu * (1 - fv) * b + (1 - fu) * fv * c + fu * fv * d
    dphi = ((1 - fv) * (b - a) + fv * (d - c)) * s
    dpsi = ((1 - fu) * (c - a) + fu * (d - b)) * s
    return lp, dphi, dpsi


def _interior_mask(n: int) -> np.ndarray:
    """Residues whose (phi, psi) PAIR both drive geometry.  `phi[0]` and `psi[n-1]` are
    never read by the builder, so index 0 and n-1 are half-inert and the pair density is
    not the right object there."""
    m = np.zeros(n, bool)
    m[1:n - 1] = True
    if not m.any():
        m[:] = True
    return m


class RamaPenalty:
    """-mean log P(phi, psi) over the constrained residues, training folds only."""

    def __init__(self, seq, fold, kind="rama"):
        self.seq, self.fold, self.kind = seq, fold, kind
        n = len(seq)
        T = logp_tables(kind)[fold]
        key = (res_classes(seq) if kind == "rama"
               else np.array([max(AA20.find(c), 0) for c in seq], int))
        self.L = T[key]
        self.mask = _interior_mask(n)
        self.nm = float(self.mask.sum())

    def __call__(self, phi, psi):
        lp = _bilinear(self.L, np.atleast_2d(phi), np.atleast_2d(psi))
        return -lp[:, self.mask].mean(1)

    def grad(self, phi, psi):
        """``(value, d/dphi, d/dpsi)`` for ONE conformation ``(n,)``."""
        lp, dp, ds = _bilinear_d(self.L, phi[None], psi[None])
        w = np.where(self.mask, -1.0 / self.nm, 0.0)
        return float(-lp[0, self.mask].mean()), w * dp[0], w * ds[0]


class RamaHingePenalty(RamaPenalty):
    """``mean_i max(0, t_class(i) - log P(phi_i, psi_i))``.

    SPARSE SUPPORT is the point.  A residue already as plausible as 95% of real residues of
    its class contributes nothing and feels no force; one in a region real protein does not
    occupy is pulled until it is merely unusual and then released.  A mean-log-density
    penalty instead has gradient everywhere and drags plausible residues toward their basin
    modes too, which is what costs it +0.008 to +0.054 A of CA-RMSD.
    """

    def __init__(self, seq, fold, kind="rama", pct=HINGE_PCT):
        super().__init__(seq, fold, kind=kind)
        self.t = hinge_thresholds(pct)[fold][res_classes(seq)]

    def __call__(self, phi, psi):
        lp = _bilinear(self.L, np.atleast_2d(phi), np.atleast_2d(psi))
        return np.maximum(self.t[None, :] - lp, 0.0)[:, self.mask].mean(1)

    def grad(self, phi, psi):
        lp, dp, ds = _bilinear_d(self.L, phi[None], psi[None])
        h = self.t - lp[0]
        act = (h > 0.0) & self.mask
        val = float(np.where(act, h, 0.0).sum() / self.nm)
        w = np.where(act, -1.0 / self.nm, 0.0)
        return val, w * dp[0], w * ds[0]


class PhiPosPenalty:
    """Smooth one-sided barrier on positive phi, non-glycine only.  The CRUDE control.

    ``sigmoid(sin(phi)/tau)`` rather than ``sigmoid(phi/tau)`` because phi lives on a
    circle: a bare logistic would read +179 deg as forbidden and -179 deg -- the same
    conformation to within 2 deg -- as allowed, and the optimiser would simply walk phi
    past +180 to escape the barrier.
    """

    def __init__(self, seq, fold, tau=PHIP_TAU):
        self.tau = tau
        m = _interior_mask(len(seq))
        mm = m & np.array([c != "G" for c in seq])
        self.mask = mm if mm.any() else m
        self.nm = float(self.mask.sum())

    def __call__(self, phi, psi):
        s = 1.0 / (1.0 + np.exp(-np.sin(np.atleast_2d(phi)) / self.tau))
        return s[:, self.mask].mean(1)

    def grad(self, phi, psi):
        sp = np.sin(phi)
        s = 1.0 / (1.0 + np.exp(-sp / self.tau))
        d = s * (1.0 - s) * np.cos(phi) / self.tau
        w = np.where(self.mask, 1.0 / self.nm, 0.0)
        return float(s[self.mask].mean()), w * d, np.zeros_like(sp)


def make_penalty(kind, seq, fold):
    """The torsion prior for one target.  Same construction as `s8.project.make_penalty`.

    `vm` is the sequence-conditioned von Mises mixture and lives in `s5.torsion` behind
    torch; it is a reporting arm, never the production penalty, so it is delegated rather
    than consolidated and it has no analytic gradient (the objective falls back to finite
    differences for it automatically).
    """
    if kind == "ramah":
        return RamaHingePenalty(seq, fold)
    if kind in ("rama", "rama20"):
        return RamaPenalty(seq, fold, kind=kind)
    if kind == "vm":
        from s8.project import VMPenalty
        return VMPenalty(seq, fold)
    if kind == "phip":
        return PhiPosPenalty(seq, fold)
    raise KeyError(kind)


# ============================================================ weighted superposition
def wkabsch_rmsd_batch(P, ref, w):
    """Weighted CA-RMSD, reflections forbidden.  At uniform weights it is numerically
    identical to the unweighted routine, which `tests/test_project.py` asserts."""
    P = np.asarray(P, float)
    ref = np.asarray(ref, float)
    if P.ndim == 2:
        P = P[None]
    w = np.asarray(w, float)
    w = w / w.sum()
    Pc = P - (w[None, :, None] * P).sum(1, keepdims=True)
    Rc = ref - (w[:, None] * ref).sum(0, keepdims=True)
    Wp = Pc * w[None, :, None]
    H = np.einsum("bni,nj->bij", Wp, Rc)
    U, S, Vt = np.linalg.svd(H)
    det = np.linalg.det(np.einsum("bji,bkj->bik", Vt, U))
    S = S.copy()
    S[:, -1] *= np.sign(det)
    resid = (Wp * Pc).sum((1, 2)) + (w[:, None] * Rc * Rc).sum() - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(resid, 0.0))


def agreement_weights(Wsub, C, floor=0.25):
    """Per-residue data weights from the local agreement of the averaged candidates."""
    d = np.sqrt(((Wsub - C[None]) ** 2).sum(-1))
    s = np.sqrt((d ** 2).mean(0))
    s0 = max(float(np.median(s)), 1e-6)
    w = 1.0 / (s + s0)
    w = w / w.mean()
    return np.maximum(w, floor)

# ============================================================ the constrained projection
#: The four generic starting conformations, in degrees: extended, alpha-helix, beta-strand,
#: polyproline II.  Copied by value from `s8.consensus2.FIT_STARTS`; the self-check asserts
#: they are still the same four.
STARTS = ((-120.0, 130.0), (-57.0, -47.0), (-139.0, 135.0), (-75.0, 145.0))


def _emit(x, n):
    """The CA trace the stage RETURNS, built by the reference builder from the final
    torsions.  The scan builder is an inner-loop device and never reaches an output."""
    geo = core.backend("geometry")
    return geo.build_backbone_batch(x[None, :n], x[None, n:])["CA"][0]


def _make_fg(C, pen, lam, w, n, grad=None):
    """The objective handed to L-BFGS-B, and whether it supplies its own gradient.

    Three modes, and which one is the default is a scientific decision, not a speed one:

    ``exact``     the reference objective, bit-for-bit: the bit-exact builder, the
                  canonical Kabsch, the reference's one-sided finite difference at the
                  reference's eps.  Reproduces `s8.project` exactly.  THE DEFAULT.
    ``analytic``  the scan builder and the true gradient.  3.3x faster again, and it does
                  NOT reproduce the reference's structures -- see `build_ca_exact`.
    ``fd``        the scan builder with the reference's finite difference.  Not for
                  production; it exists to separate "the builder moved the answer" from
                  "the gradient moved the answer", and it is what showed the builder alone
                  moves it on 126/126 targets.
    """
    mode = (grad or GRAD)
    if mode not in GRAD_MODES:
        raise KeyError(f"unknown gradient mode {mode!r}; known: {GRAD_MODES}")
    ref_c = C - C.mean(0)
    ref_sq = float((ref_c ** 2).sum())
    use_pen = bool(lam) and pen is not None
    analytic = (mode == "analytic" and w is None
                and (not use_pen or hasattr(pen, "grad")))

    if mode == "exact" and not analytic:
        #: every value on this path comes out of the same routine the reference calls
        aud = core.backend("numerics")
        E = np.eye(2 * n) * FD_EPS

        def fg(x):
            X = np.vstack([x[None], x[None] + E])
            CA = build_ca_exact(X[:, :n], X[:, n:])
            r = (aud.kabsch_rmsd_batch(CA, C) if w is None
                 else wkabsch_rmsd_batch(CA, C, w))
            if use_pen:
                r = r + lam * pen(X[:, :n], X[:, n:])
            return float(r[0]), (r[1:] - r[0]) / FD_EPS
        return fg

    if analytic:
        def fg(x):
            CA, G = build_ca(x[None, :n], x[None, n:])
            r, dca = _rmsd_and_dca(CA[0], ref_c, ref_sq)
            grad = _torsion_grad(G[0], CA[0], dca)
            if use_pen:
                pv, dp, ds = pen.grad(x[:n], x[n:])
                r = r + lam * pv
                grad[:n] += lam * dp
                grad[n:] += lam * ds
            return float(r), grad
        return fg

    #: the reference's own one-sided finite difference, over the same batched build
    E = np.eye(2 * n) * FD_EPS

    def fg(x):
        X = np.vstack([x[None], x[None] + E])
        CA = build_ca(X[:, :n], X[:, n:])[0]
        r = (_rmsd_batch(CA, ref_c, ref_sq) if w is None
             else wkabsch_rmsd_batch(CA, C, w))
        if use_pen:
            r = r + lam * pen(X[:, :n], X[:, n:])
        return float(r[0]), (r[1:] - r[0]) / FD_EPS
    return fg


def fit_prior(C, phi0, psi0, pen=None, lam=0.0, w=None, maxiter=300, grad=None):
    """The ideal-geometry chain nearest `C` that is ALSO plausible in torsion space.

        minimise   RMSD( build(phi, psi), C )  +  lam * pen(phi, psi)

    At ``lam = 0`` and ``w = None`` this is the S8-11 incumbent projection: same objective,
    same L-BFGS-B call, same cap.  Returns ``(CA, phi, psi, fval, f_start, d_to_C)``.
    """
    from scipy.optimize import minimize
    C = np.asarray(C, float)
    n = len(C)
    fg = _make_fg(C, pen, lam, w, n, grad)
    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    f0 = fg(x0)[0]
    r = minimize(fg, x0, jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
    ph, ps = r.x[:n].copy(), r.x[n:].copy()
    F = _emit(r.x, n)
    aud = core.backend("numerics")
    d = float(aud.kabsch_rmsd_batch(F[None], C)[0])
    return F, ph, ps, float(r.fun), f0, d


def fit_multi(C, pen=None, lam=0.0, w=None, extra=None, maxiter=300, grad=None):
    """`fit_prior` from every generic start plus an optional extra, keeping the best.

    The number of starts is not a tuning knob.  The projection is degenerate -- a CA trace
    admits two ideal-geometry torsion solutions at near-equal objective distance, one
    plausible and one not -- and a warm-started optimiser cannot cross between them, so a
    single start locks the answer onto whichever branch it happened to fall into.
    """
    n = len(C)
    starts = [(np.full(n, math.radians(a)), np.full(n, math.radians(b)))
              for a, b in STARTS]
    if extra is not None:
        starts.append((np.asarray(extra[0], float), np.asarray(extra[1], float)))
    best = None
    for ph, ps in starts:
        got = fit_prior(C, ph, ps, pen=pen, lam=lam, w=w, maxiter=maxiter, grad=grad)
        if best is None or got[3] < best[3]:
            best = got
    return best


def lam_path(C, pen, lams, w=None, extra=None, maxiter=300, multi=False, grad=None):
    """The path in lambda.  Solve at lam=0 from the generic starts, then walk the ladder.

    ``multi=True`` adds the generic starts back at EVERY rung and keeps the lowest
    objective, because the continuation path demonstrably under-optimises: on 1A1P, 1CS9
    and 1I6Y the multi-start finds a strictly lower objective at the same lambda with a
    positive-phi rate half as large.  The prior's real job is choosing among near-degenerate
    solutions, and it can only do that if the optimiser can see them.

    Returns ``{lam: (CA, phi, psi, fval, f_start, d_to_C)}`` with the lam=0 entry exactly
    the incumbent unconstrained `fit`.
    """
    out = {}
    cur = fit_multi(C, pen=None, lam=0.0, w=w, extra=extra, maxiter=maxiter, grad=grad)
    for lam in lams:
        if lam == 0.0:
            out[lam] = cur if w is None else fit_multi(C, pen=None, lam=0.0, w=w,
                                                       extra=extra, maxiter=maxiter,
                                                       grad=grad)
            continue
        got = fit_prior(C, cur[1], cur[2], pen=pen, lam=lam, w=w, maxiter=maxiter,
                        grad=grad)
        if multi:
            alt = fit_multi(C, pen=pen, lam=lam, w=w, extra=extra, maxiter=maxiter,
                            grad=grad)
            if alt[3] < got[3]:
                got = alt
        out[lam] = got
        cur = got
    return out

# ============================================================ self-check / bench
INPUTS_JSON = os.path.join(_ROOT, "verify", "project_inputs.json")

#: The arms `equiv` compares, in report order.  `ref` is `s8.project` itself.
ARMS = ("ref", "ex", "fd", "an")


def _outpath(name, out, limit, done=None, want=None):
    """Where a stage writes.  ONLY A complete run gets the canonical filename.

    Learned twice.  First: a one-target smoke of `equiv` silently overwrote a 126-target
    table that had cost half an hour of reference runtime, so a `limit` now goes in the
    name.  That guard was not enough, because it only covered a run LAUNCHED small.  A run
    launched full and stopped early -- by the RAM gate, which is a normal event on this box
    -- still wrote the canonical path, and a 30-target table sitting on the full table's
    filename is the worse failure of the two: it looks complete.  Aggregates read off it
    are means over an easier subset and nothing in the file says so.

    So completeness, not just intent, decides the filename.  `done < want` writes
    `..._partial{done}.json` and leaves whatever complete table exists untouched.
    """
    if out:
        return out
    stem = f"project_{name}"
    if limit:
        stem += f"_n{int(limit)}"
    elif done is not None and want is not None and done < want:
        stem += f"_partial{int(done)}"
    return os.path.join(_ROOT, "verify", stem + ".json")


def _harvested(path=None):
    p = path or INPUTS_JSON
    with open(p) as fh:
        return json.load(fh)


def harvest(manifest="tuning126", limit=None, out=None):
    """Persist the real stage-3b INPUT for every target: the coordinate average.

    Equivalence has to be measured on what the stage is actually handed, not on random
    torsions, so this runs stages 1-3a of `core.pipeline` at production settings and stores
    ``(pdb, seq, fold, n, C, native)``.  The native is a REPORTING label only -- it is
    written after the input is final and nothing in the projection ever reads it -- and it
    is here so the equivalence test can price a change in CA-RMSD rather than only in
    coordinates.
    """
    from core import pipeline as P
    cfg = P.PROD
    db = core.backend("data")
    folds = db.folds()
    rows = []
    tg = P.manifest(manifest)
    if limit:
        tg = tg[:limit]
    for t in tg:
        fold = int(folds[t.seq])
        clk = P.Clock()
        pool = P.retrieve(t, fold, cfg, clk)
        pool = P.score(pool, t.seq, fold, t.n, cfg, clk)
        got = P.filter_pool(pool, cfg, clk)
        C = P.average(got[1], got[2], clk)
        #: `filter_pool` and `average` have both grown extra return values during this
        #: sprint; take the first element when they do, so a sibling's signature change
        #: cannot silently harvest the wrong object.
        if isinstance(C, tuple):
            C = C[0]
        C = np.asarray(C, float)
        if C.shape != (int(t.n), 3):
            raise ValueError(f"average returned {C.shape}, expected {(int(t.n), 3)}")
        rows.append({"pdb": t.pdb, "seq": t.seq, "n": int(t.n), "fold": fold,
                     "C": np.asarray(C, float).tolist(),
                     "native": np.asarray(t.ca, float).tolist()})
    path = out or INPUTS_JSON
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"manifest": manifest, "config": cfg.science(), "rows": rows}, fh)
    os.replace(tmp, path)
    return rows


def geometry_report(P) -> dict:
    """Bonds, CA-CA-CA pseudo-angles, non-local clashes, CA-trace chirality, rg.

    Identical in definition to `s8.project.geometry_report` and `s8.audit8`'s, so the
    numbers here are directly comparable to the audit's table rather than merely similar.
    """
    P = np.asarray(P, float)
    n = len(P)
    step = np.linalg.norm(P[1:] - P[:-1], axis=1)
    u = P[1:] - P[:-1]
    u = u / np.linalg.norm(u, axis=1, keepdims=True)
    ang = np.degrees(np.arccos(np.clip((u[:-1] * u[1:]).sum(1), -1, 1)))
    pang = 180.0 - ang
    if n >= 4:
        b0, b1, b2 = P[1:-2] - P[:-3], P[2:-1] - P[1:-2], P[3:] - P[2:-1]
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        mm = np.cross(n1, b1 / np.linalg.norm(b1, axis=1, keepdims=True))
        dih = np.degrees(np.arctan2((mm * n2).sum(1), (n1 * n2).sum(1)))
    else:
        dih = np.array([np.nan])
    i, j = np.nonzero((np.arange(n)[None, :] - np.arange(n)[:, None]) >= 3)
    nl = np.linalg.norm(P[i] - P[j], axis=1) if len(i) else np.array([np.inf])
    d2 = ((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)
    return {"step_mean": float(step.mean()), "step_sd": float(step.std()),
            "pseudoangle_mean": float(pang.mean()),
            "frac_pseudoangle_out_of_75_150": float(((pang < 75) | (pang > 150)).mean()),
            "min_nonlocal_CA": float(nl.min()),
            "frac_nonlocal_under_4A": float((nl < 4.0).mean()),
            "frac_dihedral_positive": float(np.mean(dih > 0)),
            "rg": float(np.sqrt(d2.sum() / (2.0 * n * n)))}


def torsion_report(seq, phi, psi) -> dict:
    """Positive-phi rates.  `con` is the CONSTRAINED set: `phi[1:]`, because `phi[0]` is
    never read by the builder and reporting it measures a coin flip, not a conformation."""
    n = len(seq)
    a = (np.asarray(phi, float) + math.pi) % (2 * math.pi) - math.pi
    g = np.array([ch == "G" for ch in seq])
    con = np.zeros(n, bool)
    con[1:] = True
    ng = con & ~g
    return {"posphi_all": float((a > 0).mean()),
            "posphi_con": float((a[con] > 0).mean()) if con.any() else float("nan"),
            "posphi_con_nongly": float((a[ng] > 0).mean()) if ng.any() else float("nan"),
            "n_con_nongly": int(ng.sum()), "n_pos_con_nongly": int((a[ng] > 0).sum())}


def l_signature(phi, psi) -> float:
    """The mean sign of ``(N - CA) . [(C - CA) x (CB - CA)]`` over the emitted residues.

    This is the chirality assertion, and it is on the emitted ATOMS, not on the torsions.
    A CA trace scored by a DISTANCE is exactly mirror-blind (pinned at 1e-9 elsewhere in
    this project), so a distance objective's lowest-objective multi-start selects
    enantiomers.  This projection minimises distance to a COORDINATE cloud with reflections
    forbidden in the superposition, so the branch it picks is chirality-resolved -- but
    that safety is a property of the objective, and it would be silently lost if the
    objective were ever reformulated in distance space.

    ``+1`` IS L, calibrated on the deposited natives rather than asserted: the same triple
    product over every residue of `pdbs/1B03`, `1BTT`, `1C9A`, `1D9L`, `1DU1` and `1E0Q`
    that has all four atoms is +1 on 93 of 93 residues (mean magnitude 2.55-2.76 A^3).  A
    reflected chain returns exactly -1.
    """
    geo = core.backend("geometry")
    B = geo.build_backbone_batch(np.asarray(phi, float)[None],
                                 np.asarray(psi, float)[None])
    N, CA, C, CB = B["N"][0], B["CA"][0], B["C"][0], B["CB"][0]
    v = ((N - CA) * np.cross(C - CA, CB - CA)).sum(1)
    return float(np.sign(v).mean())


def equiv(limit=None, out=None, verbose=True):
    """Measure the optimised projection against the reference on real pipeline inputs.

    Three arms per target, all at production settings (`ramah`, lam=0.3, maxiter=300,
    multi_start=True): `s8.project` unchanged, this module's finite-difference path, and
    this module's analytic-gradient path.  Every number the sprint could have moved is
    recorded per target -- emitted coordinates, objective value, distance to the average,
    CA-RMSD to the native, torsion and geometry validity, chirality -- so a change can be
    priced instead of asserted absent.
    """
    import time
    from s8 import project as REF
    rows_in = _harvested()["rows"]
    if limit:
        rows_in = rows_in[:limit]
    aud = core.backend("numerics")

    def arms(mod, row, C, mode):
        pen = mod.make_penalty("ramah", row["seq"], int(row["fold"]))
        kw = {"grad": mode} if mode else {}
        t0 = time.perf_counter()
        p = mod.lam_path(C, pen, (0.0, 0.3), maxiter=300, multi=True, **kw)
        return p[0.0], p[0.3], time.perf_counter() - t0

    def price(tag, fit, arm, C, nat, row):
        d = {}
        d[tag + "_fval"] = float(arm[3])
        d[tag + "_fit_fval"] = float(fit[3])
        d[tag + "_dC"] = float(arm[5])
        d[tag + "_rmsd"] = float(aud.kabsch_rmsd_batch(np.asarray(arm[0])[None], nat)[0])
        d[tag + "_fit_rmsd"] = float(aud.kabsch_rmsd_batch(np.asarray(fit[0])[None], nat)[0])
        d[tag + "_L"] = l_signature(arm[1], arm[2])
        d.update({tag + "_" + k: v for k, v in geometry_report(arm[0]).items()})
        d.update({tag + "_" + k: v for k, v in torsion_report(row["seq"], arm[1],
                                                              arm[2]).items()})
        return d

    rows = []
    t_all = {k: 0.0 for k in ARMS}
    if True:
        for idx, row in enumerate(rows_in):
            C = np.asarray(row["C"], float)
            nat = np.asarray(row["native"], float)
            r = {"pdb": row["pdb"], "n": row["n"], "fold": row["fold"],
                 "avg_rmsd": float(aud.kabsch_rmsd_batch(C[None], nat)[0])}
            got = {}
            for tag, mod, mode in (("ref", REF, None),
                                   ("ex", sys.modules[__name__], "exact"),
                                   ("fd", sys.modules[__name__], "fd"),
                                   ("an", sys.modules[__name__], "analytic")):
                fit, arm, dt = arms(mod, row, C, mode)
                t_all[tag] += dt
                got[tag] = (fit, arm)
                r["t_" + tag] = dt
                r.update(price(tag, fit, arm, C, nat, row))
            for tag in ARMS[1:]:
                #: `d*` keys are the DIFFERENCE against the reference arm of the same name,
                #: kept distinct from the `*_rmsd` keys, which are distances to the native.
                for what, k in ((0, "fit"), (1, "arm")):
                    a = np.asarray(got["ref"][what][0], float)
                    b = np.asarray(got[tag][what][0], float)
                    r[f"{tag}_d{k}_maxabs"] = float(np.abs(a - b).max())
                    r[f"{tag}_d{k}_rmsd"] = float(aud.kabsch_rmsd_batch(b[None], a)[0])
                r[f"{tag}_dfval"] = r[f"{tag}_fval"] - r["ref_fval"]
                r[f"{tag}_drmsd"] = r[f"{tag}_rmsd"] - r["ref_rmsd"]
            rows.append(r)
            if verbose:
                print(f"  {idx+1:3d}/{len(rows_in)} {row['pdb']} n={row['n']:2d}  "
                      f"ref {r['t_ref']:5.2f}s  ex {r['t_ex']:5.2f}s  "
                      f"fd {r['t_fd']:5.2f}s  an {r['t_an']:5.2f}s"
                      f"   d(arm) ex {r['ex_darm_rmsd']:.1e}  fd {r['fd_darm_rmsd']:.1e}"
                      f"  an {r['an_darm_rmsd']:.1e}",
                      flush=True)
            if _mem_pct() > 92:                                     # pragma: no cover
                print("  RAM > 92%, stopping early", flush=True)
                break

    def mean(k):
        return float(np.mean([r[k] for r in rows]))

    agg = {"n_targets": len(rows), "seconds": t_all,
           "s_per_target": {k: v / max(len(rows), 1) for k, v in t_all.items()},
           "speedup": {k: t_all["ref"] / max(t_all[k], 1e-9) for k in ARMS[1:]},
           "speedup_exact": t_all["ref"] / max(t_all["ex"], 1e-9),
           "speedup_fd": t_all["ref"] / max(t_all["fd"], 1e-9),
           "speedup_analytic": t_all["ref"] / max(t_all["an"], 1e-9),
           "avg_rmsd": mean("avg_rmsd")}
    for tag in ARMS:
        agg[tag] = {"synthesis": mean(tag + "_rmsd"), "fit": mean(tag + "_fit_rmsd"),
                    "projection_cost": mean(tag + "_rmsd") - mean("avg_rmsd"),
                    "fval": mean(tag + "_fval"), "dC": mean(tag + "_dC"),
                    "step_mean": mean(tag + "_step_mean"),
                    "step_sd_max": float(max(r[tag + "_step_sd"] for r in rows)),
                    "pseudoangle_mean": mean(tag + "_pseudoangle_mean"),
                    "posphi_con_nongly": mean(tag + "_posphi_con_nongly"),
                    "clashes": float(sum(r[tag + "_frac_nonlocal_under_4A"] > 0
                                         for r in rows)),
                    "min_nonlocal_CA": float(min(r[tag + "_min_nonlocal_CA"]
                                                 for r in rows)),
                    #: +1 is L, so the WORST case is the smallest value
                    "L_worst": float(min(r[tag + "_L"] for r in rows))}
    for tag in ARMS[1:]:
        agg[tag].update({
            "darm_maxabs_worst": float(max(r[tag + "_darm_maxabs"] for r in rows)),
            "darm_rmsd_worst": float(max(r[tag + "_darm_rmsd"] for r in rows)),
            "darm_rmsd_median": float(np.median([r[tag + "_darm_rmsd"] for r in rows])),
            "dfit_rmsd_worst": float(max(r[tag + "_dfit_rmsd"] for r in rows)),
            "n_moved_1e6": int(sum(r[tag + "_darm_rmsd"] > 1e-6 for r in rows)),
            "n_moved_1e3": int(sum(r[tag + "_darm_rmsd"] > 1e-3 for r in rows)),
            "n_moved_0p01": int(sum(r[tag + "_darm_rmsd"] > 0.01 for r in rows)),
            "dfval_worst": float(max(abs(r[tag + "_dfval"]) for r in rows)),
            "dfval_mean": mean(tag + "_dfval"),
            "drmsd_mean": mean(tag + "_drmsd"),
            "drmsd_absmax": float(max(abs(r[tag + "_drmsd"]) for r in rows))})
    agg["complete"] = len(rows) == len(rows_in)
    agg["requested"] = len(rows_in)
    payload = {"agg": agg, "rows": rows}
    path = _outpath("equiv", out, limit, len(rows), len(rows_in))
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)
    return payload


def exactness(limit=None, out=None, verbose=True, equiv_path=None):
    """Does `exact` mode reproduce `s8.project` BIT-FOR-BIT on every real target?

    Checked against the reference columns already recorded in `verify/project_equiv.json`,
    which were produced by running `s8.project` itself, so no second reference run is
    needed.  The comparison is on the OBJECTIVE VALUE REACHED, the distance to the average,
    the CA-RMSD to the native and the whole geometry table -- eight scalars per target, all
    of them functions of the emitted coordinates, all required to be exactly 0.0.

    That is the end-to-end half.  The other half is `tests/test_project.py`'s check that the
    objective and its gradient are bit-identical to the reference's at arbitrary points:
    L-BFGS-B is deterministic given (f, g), so an identical (f, g) forces an identical
    trajectory and therefore an identical structure.  One is the proof, the other is the
    measurement, and neither is trusted alone.
    """
    import time
    ep = equiv_path or os.path.join(_ROOT, "verify", "project_equiv.json")
    with open(ep) as fh:
        prev = {r["pdb"]: r for r in json.load(fh)["rows"]}
    aud = core.backend("numerics")
    rows_in = _harvested()["rows"]
    if limit:
        rows_in = rows_in[:limit]
    keys = ("fval", "fit_fval", "dC", "rmsd", "fit_rmsd", "step_mean", "step_sd",
            "pseudoangle_mean", "rg", "posphi_con_nongly", "min_nonlocal_CA")
    rows, t = [], 0.0
    for idx, row in enumerate(rows_in):
        if row["pdb"] not in prev:
            continue
        C = np.asarray(row["C"], float)
        nat = np.asarray(row["native"], float)
        pen = make_penalty("ramah", row["seq"], int(row["fold"]))
        t0 = time.perf_counter()
        p = lam_path(C, pen, (0.0, 0.3), maxiter=300, multi=True, grad="exact")
        t += time.perf_counter() - t0
        fit, arm = p[0.0], p[0.3]
        got = {"fval": float(arm[3]), "fit_fval": float(fit[3]), "dC": float(arm[5]),
               "rmsd": float(aud.kabsch_rmsd_batch(np.asarray(arm[0])[None], nat)[0]),
               "fit_rmsd": float(aud.kabsch_rmsd_batch(np.asarray(fit[0])[None], nat)[0])}
        got.update(geometry_report(arm[0]))
        got.update(torsion_report(row["seq"], arm[1], arm[2]))
        d = {k: abs(got[k] - prev[row["pdb"]]["ref_" + k]) for k in keys if k in got}
        rows.append({"pdb": row["pdb"], "n": row["n"], "max_abs_diff": float(max(d.values())),
                     "per_key": {k: float(v) for k, v in d.items()},
                     "L": l_signature(arm[1], arm[2])})
        if verbose and (rows[-1]["max_abs_diff"] != 0.0 or (idx + 1) % 25 == 0):
            print(f"  {idx+1:3d}/{len(rows_in)} {row['pdb']} "
                  f"max|diff| {rows[-1]['max_abs_diff']:.3e}", flush=True)
        if _mem_pct() > 92:                                       # pragma: no cover
            break
    bad = [r for r in rows if r["max_abs_diff"] != 0.0]
    agg = {"n_targets": len(rows), "n_bit_identical": len(rows) - len(bad),
           "worst_abs_diff": float(max(r["max_abs_diff"] for r in rows)) if rows else 0.0,
           "seconds": t, "s_per_target": t / max(len(rows), 1),
           "L_worst": float(min(r["L"] for r in rows)) if rows else 0.0,
           "not_identical": [r["pdb"] for r in bad]}
    agg["complete"] = len(rows) == len(rows_in)
    agg["requested"] = len(rows_in)
    payload = {"agg": agg, "rows": rows}
    path = _outpath("exactness", out, limit, len(rows), len(rows_in))
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)
    return payload


def degeneracy(limit=None, out=None, verbose=False):
    """How big is the degeneracy?  The objective gap between the best start and the
    runner-up, per target, at the shipped rung.

    `fit_multi` returns ``argmin`` over the starts, so whenever that gap is smaller than
    the difference between two implementations of the same arithmetic, which branch is
    Returned is not determined by the objective.  This measures the gap in the objective's
    own units, so a coordinate difference between implementations can be checked against
    it instead of being called a regression or a gain.
    """
    rows_in = _harvested()["rows"]
    if limit:
        rows_in = rows_in[:limit]
    rows = []
    for row in rows_in:
        C = np.asarray(row["C"], float)
        n = len(C)
        pen = make_penalty("ramah", row["seq"], int(row["fold"]))
        vals, cas = [], []
        for aa, bb in STARTS:
            got = fit_prior(C, np.full(n, math.radians(aa)), np.full(n, math.radians(bb)),
                            pen=pen, lam=0.3, maxiter=300)
            vals.append(got[3])
            cas.append(np.asarray(got[0], float))
        o = np.argsort(vals)
        aud = core.backend("numerics")
        rows.append({"pdb": row["pdb"], "n": row["n"],
                     "best": float(vals[o[0]]),
                     "gap_to_runner_up": float(vals[o[1]] - vals[o[0]]),
                     "spread": float(np.ptp(vals)),
                     #: how far apart the two best branches are AS STRUCTURES
                     "runner_up_dist": float(aud.kabsch_rmsd_batch(
                         cas[o[1]][None], cas[o[0]])[0])})
        if verbose:
            print(f"  {row['pdb']} gap {rows[-1]['gap_to_runner_up']:.2e} "
                  f"dist {rows[-1]['runner_up_dist']:.3f}", flush=True)
    g = np.array([r["gap_to_runner_up"] for r in rows])
    dd = np.array([r["runner_up_dist"] for r in rows])
    agg = {"n_targets": len(rows), "gap_median": float(np.median(g)),
           "gap_min": float(g.min()),
           "n_gap_under_1e4": int((g < 1e-4).sum()), "n_gap_under_1e3": int((g < 1e-3).sum()),
           "n_gap_under_1e2": int((g < 1e-2).sum()),
           "runner_up_dist_median": float(np.median(dd)),
           "runner_up_dist_max": float(dd.max()),
           #: The pairing that matters: a near-tied objective whose runner-up is a
           #: DIFFERENT STRUCTURE.  A tie between two starts that converged to the same
           #: point is not a degeneracy, it is two roads to one answer -- and on this
           #: instrument every gap below 1e-4 is of that harmless kind, while the genuine
           #: branch ties sit two orders higher at 1e-3 to 1e-2.
           "n_tied_1e3_far_0p5": int(((g < 1e-3) & (dd > 0.5)).sum()),
           "n_tied_1e2_far_0p5": int(((g < 1e-2) & (dd > 0.5)).sum()),
           "n_same_point": int((dd < 1e-3).sum())}
    agg["complete"] = len(rows) == len(rows_in)
    agg["requested"] = len(rows_in)
    payload = {"agg": agg, "rows": rows}
    path = _outpath("degeneracy", out, limit, len(rows), len(rows_in))
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)
    return payload


def iteration_profile(limit=24, out=None):
    """Is `maxiter=300` binding, and does an exact gradient need fewer iterations?

    `maxiter` is a STATED PARAMETER and is not loosened anywhere here; this only measures
    whether it is reached, because "converged at 40 and kept going to 300" and "still
    moving at 300" are different stages with different levers, and nobody had looked.
    """
    from scipy.optimize import minimize as _mini
    from s8 import project as REF
    rows_in = _harvested()["rows"][:limit]
    out_rows = []
    if True:
        for mod, tag, mode in ((REF, "ref", None), (sys.modules[__name__], "ex", "exact"),
                               (sys.modules[__name__], "fd", "fd"),
                               (sys.modules[__name__], "an", "analytic")):
            kw = {"grad": mode} if mode else {}
            seen = []
            orig = _mini

            def spy(*a, **k):                                   # noqa: ANN001
                r = orig(*a, **k)
                seen.append((int(r.nit), int(r.nfev), int(r.status)))
                return r

            import scipy.optimize as so
            so.minimize = spy
            try:
                for row in rows_in:
                    pen = mod.make_penalty("ramah", row["seq"], int(row["fold"]))
                    mod.lam_path(np.asarray(row["C"], float), pen, (0.0, 0.3),
                                 maxiter=300, multi=True, **kw)
            finally:
                so.minimize = orig
            nit = np.array([s[0] for s in seen])
            nfev = np.array([s[1] for s in seen])
            out_rows.append({"arm": tag, "solves": len(seen),
                             "nit_mean": float(nit.mean()), "nit_max": int(nit.max()),
                             "nfev_mean": float(nfev.mean()), "nfev_sum": int(nfev.sum()),
                             "frac_hitting_maxiter": float((nit >= 300).mean())})
    payload = {"n_targets": len(rows_in), "arms": out_rows}
    path = _outpath("iters", out, limit)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)
    return payload


def stability(pdbs=None, limit=None, out=None, verbose=True):
    """How stable is the REFERENCE's own branch choice under an exactly null change?

    The control for the equivalence table.  Where the optimised arms land on a different
    structure it matters enormously whether that is the rewrite under-optimising or the
    stage being degenerate, and the two are told apart by asking the reference to disagree
    with ITSELF.  A CA-RMSD is exactly invariant under a rigid motion of the target cloud,
    and the penalty does not see the cloud at all, so projecting onto ``R C + t`` is
    mathematically the same problem as projecting onto ``C`` -- with every floating-point
    operation in the Kabsch differently rounded.  `s8.project` run on both is therefore the
    reference's own reproducibility, measured in the same units as the rewrite's.
    """
    import time
    from s8 import project as REF
    aud = core.backend("numerics")
    rows_in = _harvested()["rows"]
    if pdbs:
        want = set(pdbs)
        rows_in = [r for r in rows_in if r["pdb"] in want]
    if limit:
        rows_in = rows_in[:limit]
    rng = np.random.default_rng(20240904)
    rows = []
    for idx, row in enumerate(rows_in):
        C = np.asarray(row["C"], float)
        Q = np.linalg.qr(rng.normal(size=(3, 3)))[0]
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1.0                                   # a ROTATION, never a mirror
        C2 = C @ Q.T + rng.normal(scale=7.0, size=3)
        got = []
        for X in (C, C2):
            pen = REF.make_penalty("ramah", row["seq"], int(row["fold"]))
            t0 = time.perf_counter()
            p = REF.lam_path(X, pen, (0.0, 0.3), maxiter=300, multi=True)
            got.append((p[0.0], p[0.3], time.perf_counter() - t0))
        a, b = np.asarray(got[0][1][0], float), np.asarray(got[1][1][0], float)
        r = {"pdb": row["pdb"], "n": row["n"],
             "darm_rmsd": float(aud.kabsch_rmsd_batch(b[None], a)[0]),
             "dfit_rmsd": float(aud.kabsch_rmsd_batch(
                 np.asarray(got[1][0][0], float)[None], np.asarray(got[0][0][0], float))[0]),
             "dfval": float(got[1][1][3] - got[0][1][3]),
             "fval": float(got[0][1][3])}
        rows.append(r)
        if verbose:
            print(f"  {idx+1:3d}/{len(rows_in)} {row['pdb']}  self-vs-self d(arm) "
                  f"{r['darm_rmsd']:.3e}  dfval {r['dfval']:+.2e}", flush=True)
        if _mem_pct() > 92:                                       # pragma: no cover
            break
    agg = {"n_targets": len(rows),
           "darm_rmsd_worst": float(max(r["darm_rmsd"] for r in rows)),
           "darm_rmsd_median": float(np.median([r["darm_rmsd"] for r in rows])),
           "n_moved_1e3": int(sum(r["darm_rmsd"] > 1e-3 for r in rows)),
           "n_moved_0p01": int(sum(r["darm_rmsd"] > 0.01 for r in rows)),
           "n_moved_0p1": int(sum(r["darm_rmsd"] > 0.1 for r in rows)),
           "dfval_absmax": float(max(abs(r["dfval"]) for r in rows))}
    agg["complete"] = len(rows) == len(rows_in)
    agg["requested"] = len(rows_in)
    payload = {"agg": agg, "rows": rows}
    path = _outpath("stability", out, limit, len(rows), len(rows_in))
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)
    return payload


def selfcheck(verbose=True):
    """Builder against the reference, analytic gradient against central differences."""
    geo = core.backend("geometry")
    rng = np.random.default_rng(0)
    worst_build = 0.0
    for n in (3, 9, 13, 21, 40):
        phi = rng.uniform(-math.pi, math.pi, (5, n))
        psi = rng.uniform(-math.pi, math.pi, (5, n))
        ref = geo.build_backbone_batch(phi, psi)["CA"]
        got = build_ca(phi, psi)[0]
        e = float(np.abs(got - ref).max())
        worst_build = max(worst_build, e)
        if verbose:
            print(f"  build n={n:3d}  max|CA - ref| = {e:.3e} A")

    worst_grad = 0.0
    for n in (5, 12, 20):
        C = rng.normal(0.0, 5.0, (n, 3))
        seq = "".join(rng.choice(list(AA20), n))
        pen = RamaHingePenalty(seq, 0)
        x = rng.uniform(-math.pi, math.pi, 2 * n)
        ref_c = C - C.mean(0)
        ref_sq = float((ref_c ** 2).sum())
        lam = 0.3

        def f(v):
            CA = build_ca(v[None, :n], v[None, n:])[0][0]
            r = _rmsd_and_dca(CA, ref_c, ref_sq)[0]
            return r + lam * float(pen(v[:n], v[n:])[0])

        fg = _make_fg(C, pen, lam, None, n)
        ga = fg(x)[1]
        h = 1e-6
        gn = np.array([(f(x + h * np.eye(2 * n)[k]) - f(x - h * np.eye(2 * n)[k])) / (2 * h)
                       for k in range(2 * n)])
        e = float(np.abs(ga - gn).max())
        worst_grad = max(worst_grad, e)
        if verbose:
            print(f"  grad  n={n:3d}  max|analytic - central| = {e:.3e}  "
                  f"(|g| = {np.abs(gn).max():.3f})")
    if verbose:
        print(f"  STARTS {STARTS}")
    return worst_build, worst_grad


def _mem_pct():
    """RAM load via `GlobalMemoryStatusEx`.  `Get-CimInstance` costs minutes under load."""
    import ctypes

    class _MS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    m = _MS()
    m.dwLength = ctypes.sizeof(_MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))    # type: ignore[attr-defined]
    return int(m.dwMemoryLoad)


def main(argv=None):
    import time
    a = list(sys.argv[1:] if argv is None else argv)
    cmd = a[0] if a else "selfcheck"
    if cmd == "selfcheck":
        b, g = selfcheck()
        print(f"builder {b:.3e} A, gradient {g:.3e}")
        return 0
    if cmd == "harvest":
        man = a[1] if len(a) > 1 else "tuning126"
        rows = harvest(man)
        print(f"harvested {len(rows)} targets of {man} -> {INPUTS_JSON}")
        return 0
    if cmd == "equiv":
        lim = int(a[1]) if len(a) > 1 else None
        g = equiv(limit=lim)["agg"]
        print(json.dumps({k: v for k, v in g.items()}, indent=1)[:4000])
        return 0
    if cmd == "exactness":
        lim = int(a[1]) if len(a) > 1 else None
        print(json.dumps(exactness(lim)["agg"], indent=1))
        return 0
    if cmd == "degeneracy":
        lim = int(a[1]) if len(a) > 1 else None
        print(json.dumps(degeneracy(lim)["agg"], indent=1))
        return 0
    if cmd == "iters":
        lim = int(a[1]) if len(a) > 1 else 24
        print(json.dumps(iteration_profile(lim), indent=1))
        return 0
    if cmd == "stability":
        pdbs = a[1].split(",") if len(a) > 1 and not a[1].isdigit() else None
        lim = int(a[1]) if len(a) > 1 and a[1].isdigit() else None
        print(json.dumps(stability(pdbs=pdbs, limit=lim)["agg"], indent=1))
        return 0
    if cmd == "bench":
        lim = int(a[1]) if len(a) > 1 else 8
        rows = _harvested()["rows"][:lim]
        t0 = time.perf_counter()
        for r in rows:
            C = np.asarray(r["C"], float)
            pen = make_penalty("ramah", r["seq"], int(r["fold"]))
            lam_path(C, pen, (0.0, 0.3), maxiter=300, multi=True)
        dt = time.perf_counter() - t0
        print(f"GRAD={GRAD}  {dt:.3f}s / {len(rows)} targets = {dt/len(rows):.3f}s per target")
        return 0
    raise SystemExit(f"unknown command {cmd!r}")

if __name__ == "__main__":
    raise SystemExit(main())
