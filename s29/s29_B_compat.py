#!/usr/bin/env python
"""s29/s29_B_compat.py -- S29 LANE B: THE COMPATIBILITY HAMILTONIAN H = diag(E) - J M.

Pre-registered in `s29/PREREG_S29_B.md` (read it first; the falsifiers and the gate live there,
not here).  Re-opens S28-L8b / L41 under contract rule 10 with a NON-DEGENERATE off-diagonal
term, and answers the question S28 could not: with eigenvectors that carry structure, does the
state change, and does what it emits change?

THE THREE MATRICES (all native-free; all at unit spectral norm; all with zero rows/cols on the
register's padding states)
    A     the S28 Gaussian similarity graph, A_ij = exp(-d_ij^2 / 2 sigma^2), zero diagonal,
          sigma = median off-diagonal CA-RMSD (the S28 reference, not an arm).
    A_c   its double centering H A H, H = I - 11^T/N: classical MDS, so its eigenvectors are the
          pool's PRINCIPAL COORDINATES -- its principal modes of DISAGREEMENT -- and the uniform
          (typicality) mode is gone by construction, A_c 1 = 0.
    G     the SIGNED AGREEMENT matrix, G = Delta Delta^T / n_res with Delta_i member i's
          deviation from the pool mean in the readout's own frame.  G 1 = 0 because deviations
          from the mean sum to zero, and

              <psi|G|psi> = || sum_i psi_i Delta_i ||^2 / n_res,

          the squared amplitude-weighted DEPARTURE from the typical structure.  J > 0 asks for
          members whose departures reinforce, J < 0 for members whose departures cancel (which
          is what uniform averaging does).  rank(G) <= 3 n_res - 3 <= 45 EXACTLY on posed
          windows (Kabsch removes the three translations identically and the three rotations
          only to first order; S29-L11's 3 n_res - 6 <= 42 is the practical figure, not an
          identity), against the r_stable >= 266 a within-30x gradient variance would need.

THE SIGN-MIXING LEMMA (prereg section 3.2; lane T's S29-L11(d) independently): A_c 1 = G 1 = 0,
so every eigenvector with lambda != 0 is orthogonal to 1 and therefore has entries of BOTH signs.
The extremal eigenvector is a signed contrast -- one pole of a principal shape mode against the
other -- and any readout that is a function of p = psi^2 is blind to that sign and averages the
two poles back to the pool mean.  Hence the primary readout here is SIGNED (lane A's amplitude
readout, `s27/s28_A_amp.py`, reused unchanged), and the p-readout is kept to MEASURE the
cancellation rather than to hide it.

NATIVE-FREE / ORACLE SPLIT
    `gs_target_selection` sees W, E, the tie key and the matrices; it never receives a native.
    Natives are read only in `oracle_rmsd`.  `tests/test_s29_B.py` NaN-poisons `nat_ca` /
    `oracle_rr` and asserts every native-free output is bit-identical.

USAGE
    python s29/s29_B_compat.py --grad  [--limit N] [--shard i/k]   measurement 1 (independent)
    python s29/s29_B_compat.py --gs    [--limit N] [--shard i/k]   measurement 2 (ORACLE diag)
    python s29/s29_B_compat.py --analyse1                          -> s29_B_grad.json
    python s29/s29_B_compat.py --analyse2 [--targets 12|126]       -> s29_B_gs.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
GRAD_ROWS = os.path.join(RESULTS, "s29_B_grad_rows.jsonl")
GS_ROWS = os.path.join(RESULTS, "s29_B_gs_rows.jsonl")
GRAD_OUT = os.path.join(RESULTS, "s29_B_grad.json")
GS_OUT = os.path.join(RESULTS, "s29_B_gs.json")

# ---------------------------------------------------------------- settings, all pre-registered
NS = [4, 5, 6, 7, 8, 9]
N_THETA, SEED, INIT_SD, LAYERS = 120, 1009, 0.6, 3      # S28 lane B's law, so rows are comparable
ALPHA, TEMP, M_PROD = 0.18, 0.5, 75
#: prereg section 5.2 -- BOTH signs; magnitudes matched to S28's grid
J_GRID = (-3.0, -1.0, -0.3, -0.1, 0.1, 0.3, 1.0, 3.0)
#: prereg section 5.1 / 5.4 -- REAL arms, then the two controls on the two centered matrices
MAT_ARMS = ("A", "A_c", "G", "A_c|PERM", "A_c|SPEC", "G|PERM", "G|SPEC")
VAR_DIAG_N9 = 3.051e-2          # S28-L8b, `s27/results/s28_B_train.json :: full|J0` at n = 9
COH_EPS = 0.01                  # prereg A2: R4 is undefined below this sign coherence
FLOOR_A = 0.05                  # prereg A2: "the projection floor" for S29-L11 prediction 3
FAIL18 = {"1ID6", "1JBF", "1LB7", "2BFI", "2BP4", "2JN5", "2MQ2", "2N5C", "2NB7", "2NDM",
          "3BTB", "3SGO", "5W52", "7JS6", "7LCW", "8T63", "9KAR", "9L1M"}


def picks_12() -> List[str]:
    """S27 T11's 12 trainability targets -- the same set lanes T and B2 used."""
    from s25 import phys_lib as P
    return P.targets()[::11][:12]


# ============================================================ THE MATRICES (independent builders)
def unit_spectral(M: np.ndarray) -> np.ndarray:
    """M / max|lambda|. The normalisation every matrix in this lane carries."""
    M = 0.5 * (np.asarray(M, float) + np.asarray(M, float).T)
    s = float(np.abs(np.linalg.eigvalsh(M)).max())
    if s <= 0:
        raise ValueError("matrix has zero spectral norm")
    return M / s


def gaussian_similarity(D: np.ndarray, sigma: Optional[float] = None) -> Dict:
    """A_ij = exp(-d_ij^2 / 2 sigma^2), zero diagonal, sigma = median off-diagonal distance.

    Written here rather than imported so that measurement 1 is an INDEPENDENT implementation of
    lane T's rows (prereg addendum A1); `tests/test_s29_B.py` asserts it equals
    `s27.s28_B_hop.kernel_graph`'s `A` to 1e-12.
    """
    D = np.asarray(D, float)
    k = len(D)
    iu = np.triu_indices(k, 1)
    if sigma is None:
        sigma = float(np.median(D[iu]))
    A = np.exp(-(D ** 2) / (2.0 * sigma ** 2))
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return dict(A=unit_spectral(A), sigma=float(sigma))


def double_center(M: np.ndarray) -> np.ndarray:
    """H M H with H = I - 11^T/N, by row/column means (the same object, O(N^2) not O(N^3))."""
    M = np.asarray(M, float)
    r = M.mean(1, keepdims=True)
    c = M.mean(0, keepdims=True)
    return M - r - c + M.mean()


def deviations(Wp: np.ndarray) -> np.ndarray:
    """Delta: (k, 3 n_res) deviations from the pool mean of ALREADY-POSED windows `Wp`."""
    Wp = np.asarray(Wp, float)
    F = Wp.reshape(len(Wp), -1)
    return F - F.mean(0, keepdims=True)


def agreement_gram(Wp: np.ndarray) -> np.ndarray:
    """G = Delta Delta^T / n_res on posed windows. PSD, G 1 = 0, rank <= 3 n_res - 3 (exact)."""
    Delta = deviations(Wp)
    n_res = int(Wp.shape[1])
    return (Delta @ Delta.T) / float(n_res)


def perm_matrix(M: np.ndarray, perm: np.ndarray) -> np.ndarray:
    """PERM control: P M P^T. Identical spectrum; the correspondence with E is destroyed."""
    perm = np.asarray(perm, int)
    return np.asarray(M, float)[np.ix_(perm, perm)]


def spec_matrix(M: np.ndarray, rng) -> np.ndarray:
    """SPEC control (S28-L43): Q Lambda Q^T with Q Haar orthogonal and Lambda M's OWN spectrum.

    Same spectrum, no structure: the control that separates "a spread spectrum" from "these
    particular eigenvectors", which is exactly this lane's claim.
    """
    M = np.asarray(M, float)
    N = len(M)
    w = np.linalg.eigvalsh(M)
    X = rng.normal(0.0, 1.0, (N, N))
    Q, R = np.linalg.qr(X)
    Q = Q * np.sign(np.diag(R))            # the Haar fix; QR's sign convention is not uniform
    out = (Q * w) @ Q.T
    return 0.5 * (out + out.T)


def pad(M: np.ndarray, dim: int) -> np.ndarray:
    """Zero rows and columns for the register's padding states."""
    M = np.asarray(M, float)
    out = np.zeros((dim, dim))
    out[:len(M), :len(M)] = M
    return out


def spectrum_report(M: np.ndarray) -> Dict:
    """lambda_2/lambda_1, the stable rank, the top eigenvector's spread and uniform overlap."""
    M = np.asarray(M, float)
    N = len(M)
    w, V = np.linalg.eigh(M)
    o = np.argsort(-np.abs(w))
    wa, v1 = w[o], V[:, o[0]]
    fro2 = float((w ** 2).sum())
    u = np.ones(N) / np.sqrt(N)
    return dict(N=N, lam1=float(wa[0]), lam2_over_lam1=float(wa[1] / wa[0]),
                fro2=fro2, r_stable=float(fro2 / wa[0] ** 2),
                rank=int((np.abs(w) > 1e-9 * max(1.0, float(np.abs(w).max()))).sum()),
                perron_uniform_overlap2=float((v1 @ u) ** 2),
                perron_pr_over_dim=float(1.0 / (N * (v1 ** 4).sum())))


# ======================================== MEASUREMENT 1: hop-only gradient variance (independent)
def hop_grad_paramshift(circ, theta: np.ndarray, M: np.ndarray) -> np.ndarray:
    """d<psi|M|psi>/dtheta by the exact two-term shift rule. Independent of `s28_B_hop`'s."""
    S = circ.states_batch(circ._shift_grid(np.asarray(theta, float), np.pi / 2))
    vals = np.einsum("bi,ij,bj->b", S, M, S)
    return 0.5 * (vals[0::2] - vals[1::2])


def hop_grad_fd(circ, theta: np.ndarray, M: np.ndarray, h: float = 1e-5) -> np.ndarray:
    """Central finite differences: the independent check that the shift rule is the gradient."""
    S = circ.states_batch(circ._shift_grid(np.asarray(theta, float), h))
    vals = np.einsum("bi,ij,bj->b", S, M, S)
    return (vals[0::2] - vals[1::2]) / (2.0 * h)


def hop_variance(n: int, M: np.ndarray, n_theta: int = N_THETA, seed: int = SEED,
                 init_sd: float = INIT_SD, layers: int = LAYERS, fd_check: bool = False) -> Dict:
    """Var over theta ~ N(0, init_sd^2) of dF/dtheta_0 for the HOP-ONLY cost F = -<psi|M|psi>."""
    from core import quantum as Q
    circ = Q.StatevectorCircuit(n, layers)
    rng = np.random.default_rng(seed)
    P = circ.n_params()
    g0, gn, vals, fd_max = [], [], [], 0.0
    for t in range(n_theta):
        th = rng.normal(0.0, init_sd, P)
        g = -hop_grad_paramshift(circ, th, M)
        if fd_check and t < 3:
            gf = -hop_grad_fd(circ, th, M)
            den = max(float(np.abs(g).max()), 1e-12)
            fd_max = max(fd_max, float(np.abs(g - gf).max()) / den)
        g0.append(float(g[0]))
        gn.append(float(g @ g))
        vals.append(-float(circ.state(th) @ (M @ circ.state(th))))
    g0 = np.asarray(g0)
    lam1 = float(np.abs(np.linalg.eigvalsh(M)).max())
    fro2 = float((M ** 2).sum())
    return dict(n=n, layers=layers, P=P, dim=circ.dim, n_theta=n_theta,
                var_g0=float(g0.var(ddof=1)), mean_g0=float(g0.mean()),
                mean_sq_per_param=float(np.mean(gn) / P), var_F=float(np.var(vals, ddof=1)),
                lam1=lam1, fro2=fro2, r_stable=float(fro2 / max(lam1 ** 2, 1e-30)),
                pred_var=float(fro2 / circ.dim ** 2), fd_rel_max=float(fd_max))


def grad_target(pdb: str) -> List[Dict]:
    """A, A_c, G at n = 4..9 for one target. Lane T's frame convention, so the rows compare."""
    from s12 import instrument as I
    from s27 import run_pool as RP
    from s27 import s28_B_hop as B
    cand, ch, _ = RP.channels_for(pdb)
    E_full = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    order = RP.topm(E_full, cand.k, key)
    Dfull = B.pairwise_rmsd_matrix(cand.W)
    rows = []
    for n in NS:
        dim = 1 << n
        sub = order[:min(dim, cand.k)]
        Ds = Dfull[np.ix_(sub, sub)]
        g = gaussian_similarity(Ds)
        A = g["A"]
        Ac = unit_spectral(double_center(A))
        b = int(np.argmin(Ds.mean(1)))                 # lane T's medoid convention for G
        Wp = I.superpose_batch(cand.W[sub], cand.W[sub][b])
        G = unit_spectral(agreement_gram(Wp))
        for nm, M in (("A", A), ("A_c", Ac), ("G", G)):
            Mp = M if len(M) == dim else pad(M, dim)
            r = hop_variance(n, Mp, fd_check=(n == 4))
            r.update(pdb=pdb, matrix=nm, n_res=int(cand.n), k_used=int(len(sub)),
                     sigma=g["sigma"], **{("spec_" + k): v for k, v in spectrum_report(M).items()})
            rows.append(r)
    return rows


# =================================== MEASUREMENT 2: what the exact ground state selects (ORACLE)
def ground_state(E: np.ndarray, M: np.ndarray, J: float) -> Dict:
    """Exact ground state of H = diag(E) - J M by `eigh`; global sign fixed by the largest entry.

    The S28 convention fixed the sign by `sum(v) >= 0`, which is degenerate here: for a centered
    M the extremal eigenvector is orthogonal to 1, so `sum(v)` is ~0 and its sign is numerical
    noise.  The largest-|entry| convention is well defined for every state this lane produces.
    """
    Hm = np.diag(np.asarray(E, float)) - float(J) * np.asarray(M, float)
    w, V = np.linalg.eigh(Hm)
    v0 = V[:, 0]
    if v0[int(np.argmax(np.abs(v0)))] < 0:
        v0 = -v0
    p = v0 ** 2
    return dict(psi=v0, p=p / p.sum(), e0=float(w[0]), gap=float(w[1] - w[0]))


def sign_coherence(psi: np.ndarray) -> float:
    a = float(np.abs(psi).sum())
    return float((float(np.sum(psi)) / a) ** 2) if a > 0 else 0.0


def struct_diag(C: np.ndarray) -> Dict:
    C = np.asarray(C, float)
    if not np.isfinite(C).all():
        return dict(rg=float("nan"), bond=float("nan"))
    return dict(rg=float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean())),
                bond=float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean()))


def oracle_rmsd(cand, C) -> float:
    """The ONLY place a native is read. ORACLE evaluation of an ACHIEVABLE selection."""
    from s12 import instrument as I
    C = np.asarray(C, float)
    if cand.nat_ca is None or not np.isfinite(np.asarray(cand.nat_ca, float)).all():
        return float("nan")
    if not np.isfinite(C).all():
        return float("nan")
    return float(I.ca_rmsd(C, cand.nat_ca))


def gs_target_selection(W: np.ndarray, E_real: np.ndarray, key: np.ndarray, pdb: str,
                        j_grid=None, mats=None) -> Dict:
    """THE NATIVE-FREE HALF of measurement 2. Never receives a native.

    Returns {"cells": [...], "frame": Frame, "prod_idx", "C_prod", "pc1", "Delta"} -- everything
    the ORACLE half needs to score, and nothing it needs to select.
    """
    from s22 import qcand_lib as QC
    from s27 import run_pool as RP
    from s27 import s28_B_hop as B
    from s27 import s28_A_amp as AMP
    j_grid = J_GRID if j_grid is None else j_grid
    mats = MAT_ARMS if mats is None else mats
    W = np.asarray(W, float)
    E_real = np.asarray(E_real, float)
    k = len(W)
    enc = QC.Encoding(E_real)
    dim, E = enc.dim, enc.E
    top = np.sort(RP.topm(E_real, M_PROD, key)).astype(int)
    frame = AMP.Frame(W, top)                       # the deployed average's frame (top-75 medoid)
    D = B.pairwise_rmsd_matrix(W)
    Delta = deviations(frame.Wp)                    # (k, 3 n_res) in ONE frame: G, PC1, readout
    C_prod = frame.Wf[top].mean(0).reshape(frame.n, 3)
    U, S, Vt = np.linalg.svd(Delta, full_matrices=False)
    pc1 = Vt[0]                                     # (3 n_res,), unit norm: the pool's first mode
    base = {}
    A = gaussian_similarity(D)["A"]
    base["A"] = A
    base["A_c"] = unit_spectral(double_center(A))
    base["G"] = unit_spectral((Delta @ Delta.T) / float(frame.n))
    built = {}
    for name in mats:
        if "|" in name:
            root, ctrl = name.split("|")
            if ctrl == "PERM":
                built[name] = perm_matrix(base[root], RP.rng_for(pdb, "s29B_perm").permutation(k))
            elif ctrl == "SPEC":
                built[name] = spec_matrix(base[root], RP.rng_for(pdb, "s29B_spec"))
            else:
                raise ValueError(ctrl)
        else:
            built[name] = base[name]
    cells = []
    for name in mats:
        M = built[name]
        Mp = pad(M, dim)
        sr = spectrum_report(M)
        w1, V1 = np.linalg.eigh(M)
        v1 = V1[:, int(np.argmax(np.abs(w1)))]
        if v1[int(np.argmax(np.abs(v1)))] < 0:
            v1 = -v1
        for J in j_grid:
            gs = ground_state(E, Mp, J)
            psi, p = gs["psi"], gs["p"]
            psi_r, p_r = psi[:k], p[:k]
            pr_n = p_r / max(p_r.sum(), 1e-300)
            idx3 = np.sort(RP.topm(-pr_n, M_PROD, key)).astype(int)
            C3 = frame.Wf[idx3].mean(0).reshape(frame.n, 3)
            C2 = B.readout_weighted(W, D, pr_n)
            coh = sign_coherence(psi_r)
            C4, denom, wts = AMP.readout(psi_r, frame)
            if coh < COH_EPS:
                C4 = np.full((frame.n, 3), np.nan)
            eta = float((C4.reshape(-1) - C_prod.reshape(-1)) @ pc1) if np.isfinite(C4).all() else float("nan")
            if np.isfinite(eta):
                resid = C4.reshape(-1) - C_prod.reshape(-1) - eta * pc1
                num = float((C4.reshape(-1) - C_prod.reshape(-1)) @ (C4.reshape(-1) - C_prod.reshape(-1)))
                r2 = float(1.0 - (resid @ resid) / num) if num > 1e-18 else float("nan")
            else:
                resid, r2 = None, float("nan")
            sub = D[np.ix_(idx3, idx3)]
            plus = float(pr_n[v1 > 0].sum())
            minus = float(pr_n[v1 < 0].sum())
            cells.append(dict(
                matrix=name, J=float(J), e0=gs["e0"], gap=gs["gap"],
                pr=float(1.0 / np.sum(pr_n ** 2)), sign_coh=coh, denom=float(denom),
                hop=float(psi @ (Mp @ psi)),
                pole_imbalance=float(plus - minus), pole_plus=plus, pole_minus=minus,
                mean_pair_rmsd_R3=float(sub[np.triu_indices(len(idx3), 1)].mean()),
                jac75=float(len(set(idx3.tolist()) & set(top.tolist()))
                            / len(set(idx3.tolist()) | set(top.tolist()))),
                mass_top75=float(pr_n[top].sum()),
                frac_neg_w=float((wts < 0).mean()) if np.isfinite(wts).all() else float("nan"),
                ess_w=float(1.0 / max(float((wts ** 2).sum()), 1e-300)) if np.isfinite(wts).all() else float("nan"),
                eta=eta, pc1_r2=r2,
                d_R3_prod=float(np.sqrt(((C3 - C_prod) ** 2).sum(1).mean())),
                R2=dict(C=C2.tolist(), **struct_diag(C2)),
                R3=dict(idx=idx3.tolist(), C=C3.tolist(), **struct_diag(C3)),
                R4=dict(C=C4.tolist(), defined=bool(np.isfinite(C4).all()), **struct_diag(C4)),
                spec={k2: v2 for k2, v2 in sr.items()},
            ))
    return dict(cells=cells, C_prod=C_prod, pc1=pc1, top=top, frame_n=frame.n,
                base_spec={nm: spectrum_report(M) for nm, M in base.items()})


def gs_target(pdb: str) -> List[Dict]:
    """Measurement 2 for one target: the native-free selection, then the ORACLE scoring."""
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(pdb)
    E_real = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    sel = gs_target_selection(cand.W, E_real, key, pdb)
    C_prod, pc1 = sel["C_prod"], sel["pc1"]
    r_prod = oracle_rmsd(cand, C_prod)
    rows = []
    for c in sel["cells"]:
        for R in ("R2", "R3", "R4"):
            c[R]["rmsd"] = oracle_rmsd(cand, np.asarray(c[R]["C"], float))
            c[R].pop("C")
        e = c["eta"]
        if np.isfinite(e):
            plus = oracle_rmsd(cand, C_prod + abs(e) * pc1.reshape(-1, 3))
            minus = oracle_rmsd(cand, C_prod - abs(e) * pc1.reshape(-1, 3))
            c["rmsd_pc1_plus"], c["rmsd_pc1_minus"] = plus, minus
            c["rmsd_pc1_bestsign"] = float(min(plus, minus))
            c["realised_sign"] = int(np.sign(e))
            c["oracle_sign"] = int(1 if plus <= minus else -1)
            c["sign_correct"] = bool(c["realised_sign"] == c["oracle_sign"])
        c.update(pdb=pdb, n=int(cand.n), fold=int(cand.fold), k=int(cand.k),
                 basis="point_cloud", label="ACHIEVABLE_ORACLE_SCORED",
                 rmsd_prod=r_prod, fail18=bool(pdb in FAIL18),
                 prod_rg=float(np.sqrt(((C_prod - C_prod.mean(0)) ** 2).sum(1).mean())),
                 prod_bond=float(np.linalg.norm(np.diff(C_prod, axis=0), axis=1).mean()),
                 oracle_pool_best=(float(np.min(cand.oracle_rr))
                                   if cand.oracle_rr is not None
                                   and np.isfinite(np.asarray(cand.oracle_rr, float)).any()
                                   else float("nan")))
        rows.append(c)
    return rows


# ================================================================================ job plumbing
def _done(path: str) -> set:
    s = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    s.add(json.loads(line)["pdb"])
                except Exception:
                    pass
    return s


def _append(path: str, rows: List[Dict]) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    with open(path, "a", encoding="utf-8") as fh, open(tmp, encoding="utf-8") as src:
        fh.write(src.read())
    os.remove(tmp)


def _shard(pdbs: List[str], spec: Optional[str]) -> List[str]:
    if not spec:
        return pdbs
    i, k = (int(x) for x in spec.split("/"))
    return [p for j, p in enumerate(pdbs) if j % k == i]


def shard_path(path: str, spec: Optional[str]) -> str:
    """One rows file per shard: concurrent appends to one file are not atomic on Windows."""
    if not spec:
        return path
    i, k = (int(x) for x in spec.split("/"))
    return path.replace(".jsonl", f".s{i}of{k}.jsonl")


def load_all(path: str) -> List[Dict]:
    """Every shard of a rows file, de-duplicated by (pdb, and the row's own cell key)."""
    import glob
    out, seen = [], set()
    pats = [path] + sorted(glob.glob(path.replace(".jsonl", ".s*of*.jsonl")))
    for f in pats:
        if not os.path.exists(f):
            continue
        for r in load(f):
            key = (r.get("pdb"), r.get("matrix"), r.get("J"), r.get("n"),
                   r.get("lam"), r.get("kind"), r.get("greedy_m"))
            if key in seen:
                continue
            seen.add(key)
            out.append(r)
    return out


def run_phase(phase: str, pdbs: List[str], path: str, fn) -> None:
    done = _done(path)
    t0 = time.time()
    for i, pdb in enumerate(pdbs):
        if pdb in done:
            print(f"[{phase}] {pdb} done", flush=True)
            continue
        t1 = time.time()
        _append(path, fn(pdb))
        print(f"[{phase}] {i + 1}/{len(pdbs)} {pdb} {time.time() - t1:.1f}s "
              f"(elapsed {(time.time() - t0) / 60:.1f} min)", flush=True)
    print("done:", path, flush=True)


def load(path: str) -> List[Dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(l) for l in fh if l.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--grad", action="store_true")
    ap.add_argument("--gs", action="store_true")
    ap.add_argument("--analyse1", action="store_true")
    ap.add_argument("--analyse2", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--shard", type=str, default="")
    ap.add_argument("--targets", type=str, default="12")
    a = ap.parse_args(argv)
    if a.grad or a.gs:
        from s25 import phys_lib as P
        pdbs = picks_12() if a.targets == "12" else list(P.targets())
        if a.limit:
            pdbs = pdbs[:a.limit]
        pdbs = _shard(pdbs, a.shard)
        if a.grad:
            run_phase("grad", pdbs, shard_path(GRAD_ROWS, a.shard), grad_target)
        if a.gs:
            run_phase("gs", pdbs, shard_path(GS_ROWS, a.shard), gs_target)
    if a.analyse1:
        from s29 import s29_B_analyse as AN
        AN.analyse1()
    if a.analyse2:
        from s29 import s29_B_analyse as AN
        AN.analyse2(a.targets)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
