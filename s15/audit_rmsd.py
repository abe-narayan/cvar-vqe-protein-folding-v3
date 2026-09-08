"""S15 AUDIT / Part C -- an INDEPENDENT CA-RMSD implementation and its edge cases.

`s12.instrument.kabsch_rmsd_batch` computes the superposition RMSD by SVD of the
cross-covariance with the classical determinant sign fix.  The implementation here is
Horn's QUATERNION method: it builds the symmetric 4x4 key matrix K from the correlation
matrix and takes its largest eigenvalue.  That is a mathematically independent route --
different decomposition, different numerical kernel, and the rotation it implies is
PROPER BY CONSTRUCTION (a unit quaternion always maps to det = +1), so it can never be
fooled by a mirror image.  Agreement between the two is therefore evidence, not a copy.

    python -m s15.audit_rmsd            # agreement on the instrument + all edge cases
"""
from __future__ import annotations
import os, sys, json, math
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)


# ----------------------------------------------------------------- independent
def horn_rmsd(P, Q):
    """RMSD of P (n,3) onto Q (n,3) over the PROPER rotation group SO(3), Horn 1987.

    No SVD.  E0 = |P|^2 + |Q|^2 after centring; RMSD^2 = (E0 - 2 lambda_max)/n where
    lambda_max is the largest eigenvalue of the 4x4 key matrix of the correlation matrix.
    Unit quaternions cover SO(3) only, so an improper "fit" is unreachable by construction.
    """
    P = np.asarray(P, np.float64); Q = np.asarray(Q, np.float64)
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        raise ValueError("horn_rmsd wants two (n,3) arrays of equal shape")
    n = P.shape[0]
    if n == 0:
        raise ValueError("horn_rmsd: empty structure")
    P = P - P.mean(0); Q = Q - Q.mean(0)
    R = P.T @ Q                                     # 3x3 correlation matrix
    Sxx, Sxy, Sxz = R[0]; Syx, Syy, Syz = R[1]; Szx, Szy, Szz = R[2]
    K = np.array([
        [Sxx + Syy + Szz,  Syz - Szy,        Szx - Sxz,        Sxy - Syx],
        [Syz - Szy,        Sxx - Syy - Szz,  Sxy + Syx,        Szx + Sxz],
        [Szx - Sxz,        Sxy + Syx,       -Sxx + Syy - Szz,  Syz + Szy],
        [Sxy - Syx,        Szx + Sxz,        Syz + Szy,       -Sxx - Syy + Szz]],
        np.float64)
    lam = float(np.linalg.eigvalsh(K)[-1])          # symmetric -> real spectrum
    e0 = float((P * P).sum() + (Q * Q).sum())
    return math.sqrt(max(e0 - 2.0 * lam, 0.0) / n)


def horn_rotation(P, Q):
    """The optimal proper rotation (3,3) taking centred P onto centred Q, det = +1."""
    P = np.asarray(P, np.float64) - np.asarray(P, np.float64).mean(0)
    Q = np.asarray(Q, np.float64) - np.asarray(Q, np.float64).mean(0)
    R = P.T @ Q
    Sxx, Sxy, Sxz = R[0]; Syx, Syy, Syz = R[1]; Szx, Szy, Szz = R[2]
    K = np.array([
        [Sxx + Syy + Szz,  Syz - Szy,        Szx - Sxz,        Sxy - Syx],
        [Syz - Szy,        Sxx - Syy - Szz,  Sxy + Syx,        Szx + Sxz],
        [Szx - Sxz,        Sxy + Syx,       -Sxx + Syy - Szz,  Syz + Szy],
        [Sxy - Syx,        Szx + Sxz,        Syz + Szy,       -Sxx - Syy + Szz]],
        np.float64)
    w, v = np.linalg.eigh(K)
    q = v[:, -1]
    q0, q1, q2, q3 = q / np.linalg.norm(q)
    return np.array([
        [q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3),         2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3),         q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2),         2*(q2*q3+q0*q1),         q0*q0-q1*q1-q2*q2+q3*q3]])


def horn_batch(W, T):
    return np.array([horn_rmsd(w, T) for w in np.asarray(W, np.float64)])


# ---------------------------------------------------------------- edge cases
def _rand_chain(rng, n=12):
    """A random CA trace with 3.8 A virtual bonds -- physically shaped, not a Gaussian."""
    x = np.zeros((n, 3))
    d = rng.normal(size=3); d /= np.linalg.norm(d)
    for k in range(1, n):
        s = rng.normal(size=3) * 0.6
        d = d + s; d /= np.linalg.norm(d)
        x[k] = x[k - 1] + 3.8 * d
    return x


def edge_cases():
    from s12 import instrument as I
    rng = np.random.default_rng(20260905)
    out = []

    def rec(name, note, **kw):
        out.append(dict(case=name, note=note, **kw)); print(json.dumps(out[-1]), flush=True)

    # 1. identity
    A = _rand_chain(rng, 12)
    rec("identity", "same coordinates", horn=horn_rmsd(A, A), instr=I.ca_rmsd(A, A),
        expect=0.0)

    # 2. rigid motion invariance (rotation + translation)
    R = horn_rotation(_rand_chain(rng, 12), _rand_chain(rng, 12))
    B = A @ R.T + np.array([100.0, -37.5, 5.0])
    rec("rigid_motion", "rotated + translated copy must be 0",
        horn=horn_rmsd(B, A), instr=I.ca_rmsd(B, A), expect=0.0)

    # 3. MIRROR IMAGE -- the decisive test. A reflected chain is NOT superposable by a
    #    proper rotation; an implementation that forgets the determinant fix returns ~0.
    M = A.copy(); M[:, 2] *= -1.0
    naive = _naive_svd_no_sign_fix(M, A)
    rec("mirror_image", "z-reflected copy; proper-rotation RMSD must be LARGE",
        horn=horn_rmsd(M, A), instr=I.ca_rmsd(M, A), naive_svd_no_sign_fix=naive,
        expect="horn == instr >> 0, naive ~ 0")

    # 4. full inversion through the origin (det = -1, another improper map)
    rec("point_inversion", "x -> -x", horn=horn_rmsd(-A, A), instr=I.ca_rmsd(-A, A),
        naive_svd_no_sign_fix=_naive_svd_no_sign_fix(-A, A), expect="large")

    # 5. n = 1
    p = rng.normal(size=(1, 3)); q = rng.normal(size=(1, 3))
    rec("n_equals_1", "single residue: any translation fits exactly",
        horn=horn_rmsd(p, q), instr=I.ca_rmsd(p, q), expect=0.0)

    # 6. n = 2
    p = np.array([[0., 0., 0.], [3.8, 0., 0.]]); q = np.array([[1., 1., 1.], [1., 4.8, 1.]])
    rec("n_equals_2", "two residues, equal bond length: exact fit",
        horn=horn_rmsd(p, q), instr=I.ca_rmsd(p, q), expect=0.0)

    # 7. n = 3, collinear (rank-deficient correlation matrix)
    c = np.array([[0., 0., 0.], [3.8, 0., 0.], [7.6, 0., 0.]])
    c2 = c[:, [1, 0, 2]].copy()
    rec("collinear", "degenerate rank-1 geometry (rotation is not unique)",
        horn=horn_rmsd(c, c2), instr=I.ca_rmsd(c, c2), expect=0.0)

    # 8. all atoms coincident (rank-0)
    z = np.zeros((5, 3))
    rec("coincident", "every atom at the origin",
        horn=horn_rmsd(z, z + 1.0), instr=I.ca_rmsd(z, z + 1.0), expect=0.0)

    # 9. terminal residues -- full chain vs trimmed
    C = _rand_chain(rng, 14)
    D = C + rng.normal(scale=0.6, size=C.shape)
    D[0] += 6.0; D[-1] -= 6.0                       # blow up both termini only
    rec("terminal_residues", "same interior, wrecked termini: full vs [1:-1]",
        full_horn=horn_rmsd(D, C), full_instr=I.ca_rmsd(D, C),
        trimmed_horn=horn_rmsd(D[1:-1], C[1:-1]), trimmed_instr=I.ca_rmsd(D[1:-1], C[1:-1]),
        expect="trimmed MUCH lower -- this is the definition the preprint left unstated")

    # 10. missing atoms / NaN propagation
    E = C.copy(); E[3] = np.nan
    try:
        h = horn_rmsd(E, C)
    except Exception as ex:
        h = f"{type(ex).__name__}"
    try:
        iv = I.ca_rmsd(E, C)
    except Exception as ex:
        iv = type(ex).__name__
    rec("missing_atom_nan", "a NaN coordinate must NOT silently produce a finite number",
        horn=h, instr=iv, expect="nan or raise, never a plausible float")

    # 11. mismatched lengths
    try:
        h = horn_rmsd(C[:10], C)
    except Exception as ex:
        h = type(ex).__name__
    try:
        v = I.ca_rmsd(C[:10], C); v = float(v)
    except Exception as ex:
        v = type(ex).__name__
    rec("length_mismatch", "10 vs 14 residues", horn=h, instr=v, expect="raise")

    # 12. chain break -- geometry the builder can never emit
    F = C.copy(); F[7:] += 25.0
    rec("chain_break", "a 25 A discontinuity mid-chain is scored, not rejected",
        horn=horn_rmsd(F, C), instr=I.ca_rmsd(F, C), expect="finite, large")

    # 13. units / scale sanity: doubling every coordinate
    rec("scale_2x", "RMSD is in the coordinate unit (Angstrom); it is NOT scale invariant",
        horn=horn_rmsd(2 * C, C), instr=I.ca_rmsd(2 * C, C))

    # 14. permutation (alternate conformation / atom-order swap)
    G = C.copy(); G[[4, 5]] = G[[5, 4]]
    rec("atom_order_swap", "swapping two residues changes the value: correspondence is by "
        "INDEX, never by nearest neighbour", horn=horn_rmsd(G, C), instr=I.ca_rmsd(G, C))

    # 15. float32 vs float64 input
    rec("float32_input", "the universes store W as float32",
        horn64=horn_rmsd(C, D), horn32=horn_rmsd(C.astype(np.float32).astype(np.float64),
                                                 D.astype(np.float32).astype(np.float64)))
    return out


def _naive_svd_no_sign_fix(P, Q):
    """The WRONG implementation: SVD superposition without the determinant correction.
    Included only to show what the mirror test is meant to catch."""
    P = np.asarray(P, float) - np.asarray(P, float).mean(0)
    Q = np.asarray(Q, float) - np.asarray(Q, float).mean(0)
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    num = (P ** 2).sum() + (Q ** 2).sum() - 2.0 * S.sum()
    return math.sqrt(max(num, 0.0) / len(P))


# ------------------------------------------------------------ instrument check
def agreement(nmax=126):
    """Horn vs the instrument on real pool geometry: every target, K=500 windows."""
    from s12 import instrument as I
    worst = 0.0; worst_pdb = None; rows = []
    for t in I.targets()[:nmax]:
        u = I.load_univ(t["pdb"])
        p = I.pool_idx(u)
        W = u["W"][p]
        a = I.kabsch_rmsd_batch(W, u["nat_ca"])
        b = horn_batch(W, u["nat_ca"])
        d = float(np.abs(a - b).max())
        # and against the STORED oracle labels rr
        rrd = float(np.abs(np.asarray(u["rr"], np.float64)[p] - b).max())
        rows.append({"pdb": t["pdb"], "n": t["n"], "max_abs_diff_horn_vs_instr": d,
                     "max_abs_diff_horn_vs_stored_rr": rrd})
        if d > worst:
            worst, worst_pdb = d, t["pdb"]
    res = {"n_targets": len(rows), "n_structures": 500 * len(rows),
           "max_abs_diff_horn_vs_instrument": worst, "at": worst_pdb,
           "max_abs_diff_horn_vs_stored_rr": max(r["max_abs_diff_horn_vs_stored_rr"] for r in rows),
           "note_stored_rr": "s8/generate_univ stores rr as float32; ~1e-7 relative is expected",
           "rows": rows}
    return res


if __name__ == "__main__":
    res = {"edge_cases": edge_cases()}
    print("--- agreement on the instrument ---", flush=True)
    res["agreement"] = agreement()
    a = res["agreement"]
    print(json.dumps({k: a[k] for k in a if k != "rows"}, indent=1))
    with open(os.path.join(OUT, "audit_rmsd.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print("wrote", os.path.join(OUT, "audit_rmsd.json"))
