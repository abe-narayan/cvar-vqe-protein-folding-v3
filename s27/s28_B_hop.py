#!/usr/bin/env python
"""s27/s28_B_hop.py -- S28 LANE B: THE FIRST NON-DIAGONAL HAMILTONIAN, H = diag(E) - J A.

Pre-registered in `s27/PREREG_S28_B.md` (read it first; the falsifiers live there, not here).

WHAT IS NEW
===========
Every Hamiltonian this project has run is diagonal in the candidate basis, and the S24/S25
set-equality theorem makes the CVaR tail of any such H a subset of the classical energy prefix.
Here A is the pool's structural-similarity graph (Gaussian kernel on the pairwise CA-RMSD of the
500 windows, native-free) and H = diag(E) - J A is a tight-binding Hamiltonian on the pool. The
objective keeps the spine,

    F(theta) = CVaR_alpha(E; p_theta) - T H(p_theta) - J <psi_theta|A|psi_theta>,

with the first two terms and their gradient taken from `core.quantum.free_energy` UNCHANGED and
the hopping term evaluated on the real statevector with its own exact parameter-shift gradient.
At J = 0 the training loop is `core.quantum.run_cvar_vqe` line for line (same RNG, same Adam),
so J = 0 reproduces `s24.d_harness.arm_vqe` bit for bit (tested).

For J > 0 the objective is NOT a function of p alone: A >= 0 entrywise, so two states with the
same p differ in F through the signs of their amplitudes. Every earlier objective here was a
function of p.

NATIVE-FREE / ORACLE SPLIT
==========================
`select_target` sees `W`, `E`, the tie key and the graphs; it never receives a native. Natives
are read only in `oracle_rmsd`. `tests/test_s28_B.py` NaN-poisons `nat_ca` / `oracle_rr` and
asserts the selection output is bit-identical.

USAGE
=====
    python s27/s28_B_hop.py --run [--limit N]         point-cloud endpoint, all arms, resumable
    python s27/s28_B_hop.py --train                    F5: gradient variance vs J at n = 4..9
    python s27/s28_B_hop.py --chain --arms SPEC        built chain for named arms
    python s27/s28_B_hop.py --analyse                  statistics -> s28_B_summary.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q              # noqa: E402
from s12 import instrument as I            # noqa: E402
from s22 import qcand_lib as QC            # noqa: E402
from s24 import d_harness as H             # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s28_B_rows.jsonl")
CHAIN_ROWS = os.path.join(RESULTS, "s28_B_chain_rows.jsonl")
TRAIN_OUT = os.path.join(RESULTS, "s28_B_train.json")
SUMMARY = os.path.join(RESULTS, "s28_B_summary.json")

#: S25 / S27 selector settings, held fixed (`s27/PREREG_S28_B.md` section 1).
ALPHA, TEMP, LAYERS, ITERS, LR = 0.18, 0.5, 3, 80, 0.15
M_PROD = 75
J_GRID = (0.0, 0.1, 0.3, 1.0, 3.0)
GRAPHS = ("REAL", "PERM", "RAND")
SEEDS = (0, 1)
FAIL18 = {"1ID6", "1JBF", "1LB7", "2BFI", "2BP4", "2JN5", "2MQ2", "2N5C", "2NB7", "2NDM",
          "3BTB", "3SGO", "5W52", "7JS6", "7LCW", "8T63", "9KAR", "9L1M"}


# ================================================================= the pool graph (native-free)
def pairwise_rmsd_matrix(W: np.ndarray) -> np.ndarray:
    """(k, k) CA-RMSD after optimal superposition, symmetrised, zero diagonal.

    `core.backend("numerics").kabsch_rmsd_batch` gives one row per call (the brief's route);
    the matrix is symmetrised because the two one-sided superpositions agree only to float
    round-off.
    """
    import core
    kb = core.backend("numerics").kabsch_rmsd_batch
    W = np.asarray(W, float)
    k = len(W)
    D = np.empty((k, k))
    for a in range(k):
        D[a] = kb(W, W[a])
    D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    return D


def kernel_graph(D: np.ndarray, sigma: Optional[float] = None) -> Dict:
    """A_raw_ij = exp(-d_ij^2 / 2 sigma^2), zero diagonal; sigma = median off-diagonal d.

    Returns the raw graph, its spectral radius, the unit-spectral-norm graph `A` (the one the
    Hamiltonian uses), the degree sequence and sigma.
    """
    D = np.asarray(D, float)
    k = len(D)
    iu = np.triu_indices(k, 1)
    if sigma is None:
        sigma = float(np.median(D[iu]))
    A_raw = np.exp(-D ** 2 / (2.0 * sigma ** 2))
    np.fill_diagonal(A_raw, 0.0)
    return finish_graph(A_raw, sigma)


def finish_graph(A_raw: np.ndarray, sigma: float) -> Dict:
    A_raw = 0.5 * (A_raw + A_raw.T)
    np.fill_diagonal(A_raw, 0.0)
    lam = np.linalg.eigvalsh(A_raw)
    lam_max = float(lam[-1])
    return dict(A_raw=A_raw, A=A_raw / lam_max, lam_max=lam_max, lam_2=float(lam[-2]),
                lam_min=float(lam[0]), degree=A_raw.sum(1), sigma=float(sigma),
                mean_offdiag=float(A_raw[np.triu_indices(len(A_raw), 1)].mean()))


def permuted_graph(g: Dict, perm: np.ndarray) -> Dict:
    """PERM control: P A P^T. Same spectrum and degree multiset; correspondence with E destroyed."""
    perm = np.asarray(perm, int)
    A = g["A_raw"][np.ix_(perm, perm)]
    out = finish_graph(A, g["sigma"])
    out["perm"] = perm
    return out


def random_degree_matched_graph(g: Dict, rng: np.random.Generator, dim_embed: int = 3,
                                tol: float = 1e-8, max_iter: int = 20000) -> Dict:
    """RAND control: a random geometric graph with the SAME DEGREE SEQUENCE as `g`.

    500 points iid N(0, I_3), the same Gaussian kernel with its own median-distance sigma, then
    symmetric Sinkhorn scaling x_i R_ij x_j so that every row sum equals the real graph's row
    sum to `tol` (relative to the largest degree). Candidate i keeps its own degree; only who is
    similar to whom is randomised.
    """
    d = np.asarray(g["degree"], float)
    k = len(d)
    X = rng.normal(0.0, 1.0, (k, dim_embed))
    Dr = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    iu = np.triu_indices(k, 1)
    s = float(np.median(Dr[iu]))
    R = np.exp(-Dr ** 2 / (2.0 * s ** 2))
    np.fill_diagonal(R, 0.0)
    x = np.ones(k)
    scale = float(d.max())
    ok = False
    for it in range(max_iter):
        y = R @ x
        err = float(np.max(np.abs(x * y - d))) / scale
        if err < tol:
            ok = True
            break
        x = np.sqrt(x * d / np.maximum(y, 1e-300))
    A = (x[:, None] * R) * x[None, :]
    np.fill_diagonal(A, 0.0)
    out = finish_graph(A, s)
    out["sinkhorn_converged"] = bool(ok)
    out["sinkhorn_iters"] = int(it + 1)
    out["sinkhorn_err"] = float(err)
    out["max_entry"] = float(A.max())
    return out


def pad_graph(A: np.ndarray, dim: int) -> np.ndarray:
    """Zero rows and columns for the padding states of the register."""
    k = len(A)
    Ap = np.zeros((dim, dim))
    Ap[:k, :k] = A
    return Ap


# ============================================================ the hopping term and its gradient
def hop_value(psi: np.ndarray, A: np.ndarray) -> float:
    """<psi|A|psi> on a real statevector."""
    return float(psi @ (A @ psi))


def hop_abs_bound(psi: np.ndarray, A: np.ndarray) -> float:
    """sum_ij |psi_i| A_ij |psi_j|: the value the same p would reach with aligned signs."""
    a = np.abs(psi)
    return float(a @ (A @ a))


def sign_coherence(psi: np.ndarray) -> float:
    """(sum psi)^2 / (sum |psi|)^2 in [0, 1]; 1 when every amplitude has the same sign."""
    s = float(np.abs(psi).sum())
    return float((psi.sum() / s) ** 2) if s > 0 else 0.0


def grad_hop_paramshift(circ: Q.StatevectorCircuit, theta: np.ndarray, A: np.ndarray
                        ) -> np.ndarray:
    """d <psi|A|psi> / d theta_k by the exact two-term shift rule (pi/2, factor 1/2).

    `<psi|A|psi>` is the expectation of a Hermitian observable in a state prepared by RY and
    CNOT gates, so each parameter obeys the standard shift rule. All 2P shifted STATES (not
    probabilities) are simulated in one batched pass.
    """
    theta = np.asarray(theta, float)
    S = circ.states_batch(circ._shift_grid(theta, np.pi / 2))      # (2P, dim)
    AS = S @ A                                                     # A symmetric: (2P, dim)
    vals = (AS * S).sum(1)
    return (vals[0::2] - vals[1::2]) / 2.0


def grad_hop_fd(circ: Q.StatevectorCircuit, theta: np.ndarray, A: np.ndarray,
                h: float = 1e-5) -> np.ndarray:
    """Central finite differences on <psi|A|psi>: the independent check on the shift rule."""
    theta = np.asarray(theta, float)
    S = circ.states_batch(circ._shift_grid(theta, h))
    vals = np.einsum("bi,ij,bj->b", S, A, S)
    return (vals[0::2] - vals[1::2]) / (2 * h)


def hop_objective(circ: Q.StatevectorCircuit, theta: np.ndarray, E: np.ndarray, A: np.ndarray,
                  alpha: float, T: float, J: float):
    """F = CVaR_alpha(E; p) - T H(p) - J <psi|A|psi> and its exact gradient.

    Returns (F, grad F, p, CVaR, H, hop). At J == 0 the gradient is `free_energy`'s, untouched
    (no zero-times-array is added, so it is bit-identical).
    """
    f, g, p, cv, Hent = Q.free_energy(circ, theta, E, alpha, T)
    if J == 0.0:
        return float(f), g, p, float(cv), float(Hent), 0.0
    psi = circ.state(theta)
    hv = hop_value(psi, A)
    gh = grad_hop_paramshift(circ, theta, A)
    return float(f - J * hv), g - J * gh, p, float(cv), float(Hent), float(hv)


def run_hop_vqe(E: np.ndarray, A: np.ndarray, J: float, alpha: float = ALPHA, T: float = TEMP,
                n: int = 9, layers: int = LAYERS, iters: int = ITERS, restarts: int = 1,
                seed: int = 0, lr: float = LR):
    """`core.quantum.run_cvar_vqe` with the hopping term. Same RNG draw, same Adam, same order
    of operations; at J = 0 every float is identical to the deployed loop.

    Returns (p, cvar, H, hop, F, theta, circ).
    """
    circ = Q.StatevectorCircuit(n, layers)
    rng = np.random.default_rng(seed)
    best, best_f = None, np.inf
    for r in range(restarts):
        th = rng.normal(0.0, 0.6, circ.n_params())
        m = np.zeros_like(th)
        v = np.zeros_like(th)
        for t in range(1, iters + 1):
            f, g, _, _, _, _ = hop_objective(circ, th, E, A, alpha, T, J)
            m = 0.9 * m + 0.1 * g
            v = 0.999 * v + 0.001 * g * g
            th = th - lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
        f, _, p, cv, Hent, hv = hop_objective(circ, th, E, A, alpha, T, J)
        if f < best_f:
            best_f, best = f, (p, cv, Hent, hv, th.copy())
    p, cv, Hent, hv, th = best
    return p, float(cv), float(Hent), float(hv), float(best_f), th, circ


# ============================================================================ exact ground state
def ground_state(E: np.ndarray, A: np.ndarray, J: float) -> Dict:
    """The exact ground state of H = diag(E) - J A by `numpy.linalg.eigh`; p = v^2.

    At J = 0 H is diagonal and the ground state is the one-hot argmin (degenerate readouts,
    labelled so by the caller).
    """
    Hm = np.diag(np.asarray(E, float)) - float(J) * np.asarray(A, float)
    w, v = np.linalg.eigh(Hm)
    v0 = v[:, 0]
    if v0.sum() < 0:                      # fix the global sign (Perron-like states are positive)
        v0 = -v0
    p = v0 ** 2
    p = p / p.sum()
    return dict(p=p, psi=v0, e0=float(w[0]), gap=float(w[1] - w[0]))


# ================================================================================= readouts
def readout_tail(enc: QC.Encoding, p: np.ndarray, alpha: float) -> Dict:
    """R1: the alpha-tail support of p along the E order (the deployed readout), real
    candidates. `exact_face` is the S22 object; nothing here is re-implemented."""
    face = QC.exact_face(enc.E, np.asarray(p, float), alpha)
    cands = QC.tail_candidates(enc, face["tail"])
    return dict(idx=np.sort(cands).astype(int), m=int(len(cands)),
                tail_size=int(face["tail_size"]), ess=float(face["ess"]),
                n_pad=int(face["tail_size"] - len(cands)))


def readout_weighted(W: np.ndarray, D: np.ndarray, w: np.ndarray) -> np.ndarray:
    """R2: the p-weighted coordinate average in the p-weighted consensus-medoid frame
    (`core.pipeline.consensus_medoid`, the S23 L8 / `core.pipeline.average_weighted`
    convention) over ALL real candidates."""
    from core import pipeline as PL
    w = np.asarray(w, float)
    w = w / w.sum()
    b = PL.consensus_medoid(D, w)
    Sup = I.superpose_batch(np.asarray(W, float), np.asarray(W, float)[b])
    return np.tensordot(w, Sup, axes=(0, 0))


def readout_ptop(p_real: np.ndarray, m: int, key: np.ndarray) -> np.ndarray:
    """R3: the m most probable real candidates, exact ties broken by the stable key
    (`s27.run_pool.topm`), never by array order."""
    from s27 import run_pool as RP
    return np.sort(RP.topm(-np.asarray(p_real, float), int(m), key)).astype(int)


def set_descriptors(idx: np.ndarray, D: np.ndarray, dis_top: np.ndarray) -> Dict:
    """Native-free descriptors of a retained set: its mean pairwise RMSD (the cluster it
    selects) and its overlap with the DIS top-75."""
    idx = np.asarray(idx, int)
    if len(idx) > 1:
        sub = D[np.ix_(idx, idx)]
        mpr = float(sub[np.triu_indices(len(idx), 1)].mean())
    else:
        mpr = 0.0
    ov = float(len(set(idx.tolist()) & set(np.asarray(dis_top, int).tolist())) / max(1, len(dis_top)))
    return dict(mean_pair_rmsd=mpr, overlap_dis75=ov)


# ===================================================================== per-target, native-free
def build_graphs(W: np.ndarray, pdb: str) -> Tuple[np.ndarray, Dict[str, Dict]]:
    """The pairwise matrix and the three graphs (REAL, PERM, RAND) for one pool."""
    from s27 import run_pool as RP
    D = pairwise_rmsd_matrix(W)
    real = kernel_graph(D)
    k = len(W)
    perm = RP.rng_for(pdb, "s28B_perm").permutation(k)
    graphs = {"REAL": real, "PERM": permuted_graph(real, perm),
              "RAND": random_degree_matched_graph(real, RP.rng_for(pdb, "s28B_rand"))}
    return D, graphs


def consume(p: np.ndarray, psi: Optional[np.ndarray], enc: QC.Encoding, W: np.ndarray,
            D: np.ndarray, key: np.ndarray, dis_top: np.ndarray, A_pad: np.ndarray,
            p_ref: Optional[np.ndarray]) -> Dict:
    """All three readouts and the departure diagnostics of one probability vector. Native-free."""
    k = len(W)
    p = np.asarray(p, float)
    p_real = p[:k]
    p_real_n = p_real / p_real.sum()
    r1 = readout_tail(enc, p, ALPHA)
    C1, _ = I.coordinate_average(W[r1["idx"]])      # == s24.d_harness.readout_uniform's operator
    C2 = readout_weighted(W, D, p_real_n)
    r3 = readout_ptop(p_real_n, M_PROD, key)
    C3, _ = I.coordinate_average(W[r3])
    gate = H.gate_set_equality(enc.E[:k], r1["idx"], r1["m"])
    dis_set = set(np.asarray(dis_top, int).tolist())
    r3_set = set(r3.tolist())
    out = dict(
        m=r1["m"], tail_size=r1["tail_size"], n_pad_in_tail=r1["n_pad"], ess=r1["ess"],
        gate_pass=bool(gate["pass_"]), gate_equality=bool(gate["equality"]),
        gate_n_holes=int(gate["n_holes"]),
        jac75=float(len(r3_set & dis_set) / len(r3_set | dis_set)),
        mass_top75=float(p_real_n[list(dis_set)].sum()),
        pad_mass=float(p[k:].sum()),
        pr=float(1.0 / np.sum(p_real_n ** 2)),
        entropy_bits=float(-(p[p > 0] * np.log2(p[p > 0])).sum()),
        tv_vs_J0=float(0.5 * np.abs(p - p_ref).sum()) if p_ref is not None else None,
        R1=dict(idx=r1["idx"].tolist(), C=C1.tolist(), **set_descriptors(r1["idx"], D, dis_top)),
        R2=dict(C=C2.tolist()),
        R3=dict(idx=r3.tolist(), C=C3.tolist(), **set_descriptors(r3, D, dis_top)),
    )
    if psi is not None:
        out["hop"] = hop_value(psi, A_pad)
        out["hop_abs"] = hop_abs_bound(psi, A_pad)
        out["sign_coh"] = sign_coherence(psi)
    return out


def select_target(W: np.ndarray, E: np.ndarray, key: np.ndarray, pdb: str,
                  j_grid=None, graphs=None, seeds=None, with_gs: bool = True) -> List[Dict]:
    """THE NATIVE-FREE HALF. Every arm's retained sets, coordinates and diagnostics for one
    pool. Never receives a native. Returns one row per (source, seed, graph, J)."""
    from s27 import run_pool as RP
    j_grid = J_GRID if j_grid is None else j_grid       # resolved at call time (tests patch them)
    graphs = GRAPHS if graphs is None else graphs
    seeds = SEEDS if seeds is None else seeds
    W = np.asarray(W, float)
    E = np.asarray(E, float)
    k = len(W)
    enc = QC.Encoding(E)                              # identity label: candidate i at bit i
    dim = enc.dim
    dis_top = RP.topm(E, M_PROD, key)
    D, G = build_graphs(W, pdb)
    rows = []
    p_ref = {}
    # J = 0 once per source (the graph is multiplied by zero); graph tag NONE
    zero = np.zeros((dim, dim))
    for s in seeds:
        t0 = time.time()
        p, cv, Hent, hv, F, th, circ = run_hop_vqe(enc.E, zero, 0.0, seed=s, n=enc.n_qubits)
        psi = circ.state(th)
        p_ref[("vqe", s)] = p
        row = dict(source="vqe", seed=s, graph="NONE", J=0.0, cvar=cv, entropy_nats=Hent, F=F,
                   secs=time.time() - t0)
        row.update(consume(p, psi, enc, W, D, key, dis_top, zero, None))
        rows.append(row)
    if with_gs:
        gs = ground_state(enc.E, zero, 0.0)
        p_ref[("gs", -1)] = gs["p"]
        row = dict(source="gs", seed=-1, graph="NONE", J=0.0, e0=gs["e0"], gap=gs["gap"],
                   degenerate=True, secs=0.0)
        row.update(consume(gs["p"], gs["psi"], enc, W, D, key, dis_top, zero, None))
        rows.append(row)
    for gname in graphs:
        g = G[gname]
        A_pad = pad_graph(g["A"], dim)
        ginfo = dict(sigma=g["sigma"], lam_max=g["lam_max"], lam_2=g["lam_2"],
                     lam_min=g["lam_min"], mean_offdiag=g["mean_offdiag"],
                     deg_E_corr=float(np.corrcoef(g["degree"], E)[0, 1]))
        if gname == "RAND":
            ginfo.update(sinkhorn_converged=g["sinkhorn_converged"],
                         sinkhorn_err=g["sinkhorn_err"], rand_max_entry=g["max_entry"])
        for J in j_grid:
            if J == 0.0:
                continue
            for s in seeds:
                t0 = time.time()
                p, cv, Hent, hv, F, th, circ = run_hop_vqe(enc.E, A_pad, J, seed=s,
                                                           n=enc.n_qubits)
                psi = circ.state(th)
                row = dict(source="vqe", seed=s, graph=gname, J=float(J), cvar=cv,
                           entropy_nats=Hent, F=F, secs=time.time() - t0, **ginfo)
                row.update(consume(p, psi, enc, W, D, key, dis_top, A_pad, p_ref[("vqe", s)]))
                rows.append(row)
            if with_gs:
                t0 = time.time()
                gs = ground_state(enc.E, A_pad, J)
                row = dict(source="gs", seed=-1, graph=gname, J=float(J), e0=gs["e0"],
                           gap=gs["gap"], degenerate=False, secs=time.time() - t0, **ginfo)
                row.update(consume(gs["p"], gs["psi"], enc, W, D, key, dis_top, A_pad,
                                   p_ref[("gs", -1)]))
                rows.append(row)
    return rows


# ======================================================================== ORACLE scoring
def oracle_rmsd(cand, C) -> float:
    """The ONLY place a native is read. ORACLE evaluation of an ACHIEVABLE selection."""
    if cand.nat_ca is None or not np.isfinite(np.asarray(cand.nat_ca, float)).all():
        return float("nan")                       # a NaN-poisoned native scores NaN, never raises
    return float(I.ca_rmsd(np.asarray(C, float), cand.nat_ca))


def run_target(pdb: str, with_gs: bool = True) -> List[Dict]:
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(pdb)
    E = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    rows = select_target(cand.W, E, key, pdb, with_gs=with_gs)
    dis_top = RP.topm(E, M_PROD, key)
    r_dis75 = oracle_rmsd(cand, I.coordinate_average(cand.W[dis_top])[0])
    for r in rows:
        r.update(pdb=pdb, n=int(cand.n), fold=int(cand.fold), k=int(cand.k),
                 basis="point_cloud", label="ACHIEVABLE", rmsd_dis75=r_dis75,
                 fail18=bool(pdb in FAIL18))
        for R in ("R1", "R2", "R3"):
            r[R]["rmsd"] = oracle_rmsd(cand, np.asarray(r[R]["C"], float))
        r["oracle_pool_best"] = float(np.min(cand.oracle_rr)) if cand.oracle_rr is not None else float("nan")
    return rows


def arm_of(r: Dict) -> str:
    return f"{r['source']}|s{r['seed']}|{r['graph']}|J{r['J']:g}"


# ========================================================== F5: trainability vs J (no RMSD)
def deployed_E(dim: int) -> np.ndarray:
    """`s25/q_plateau.py :: deployed_E`: the standardised rank ladder over `dim` states."""
    r = np.arange(1, dim + 1, dtype=float)
    return (r - r.mean()) / r.std()


def measure_hop(n: int, layers: int, alpha: float, T: float, J: float, n_theta: int, seed: int,
                init_sd: float, E: np.ndarray, A: np.ndarray, mode: str = "full") -> Dict:
    """`s25/q_plateau.py :: measure` extended with the hopping term.

    mode = "full"      F = CVaR_alpha(E; p) - T H(p) - J <psi|A|psi>
    mode = "hop_only"  F = -<psi|A|psi>            (the linear cost with a NON-diagonal observable)
    mode = "linear"    F = <psi|diag(E)|psi>       (alpha = 1, T = 0, J = 0: the S25 control)
    """
    circ = Q.StatevectorCircuit(n, layers)
    rng = np.random.default_rng(seed)
    P = circ.n_params()
    g0, gn, vals = [], [], []
    for _ in range(n_theta):
        th = rng.normal(0.0, init_sd, P)
        if mode == "hop_only":
            g = -grad_hop_paramshift(circ, th, A)
            f = -hop_value(circ.state(th), A)
        elif mode == "linear":
            g = Q.grad_cvar_paramshift(circ, th, E, 1.0)
            f = float(circ.probs(th) @ E)
        else:
            f, g, _, _, _, _ = hop_objective(circ, th, E, A, alpha, T, J)
        g0.append(float(g[0]))
        gn.append(float(np.dot(g, g)))
        vals.append(float(f))
    g0 = np.asarray(g0)
    gn = np.asarray(gn)
    return dict(n=n, layers=layers, P=P, dim=circ.dim, alpha=alpha, T=T, J=float(J), mode=mode,
                var_g0=float(g0.var(ddof=1)), mean_g0=float(g0.mean()),
                mean_sq_norm=float(gn.mean()), mean_sq_per_param=float(gn.mean() / P),
                var_F=float(np.var(vals, ddof=1)), n_theta=n_theta)


def train_main(n_theta: int = 120):
    from s25 import phys_lib as P
    from s27 import run_pool as RP
    pdbs = P.targets()
    pick = pdbs[::11][:12]                                # S27 T11's 12 targets
    ns = [4, 5, 6, 7, 8, 9]
    rows = []
    t0 = time.time()
    for pdb in pick:
        cand, ch, _ = RP.channels_for(pdb)
        E_full = RP.zr(ch["DIS"])
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        order = RP.topm(E_full, cand.k, key)
        D = pairwise_rmsd_matrix(cand.W)
        for n in ns:
            dim = 1 << n
            if dim <= cand.k:
                sub = order[:dim]
                g = kernel_graph(D[np.ix_(sub, sub)])
                A = g["A"]
                E = deployed_E(dim)
            else:                                          # n = 9: 500 real + 12 padding
                g = kernel_graph(D)
                A = pad_graph(g["A"], dim)
                enc = QC.Encoding(E_full)
                E = enc.E
            for J in J_GRID:
                r = measure_hop(n, LAYERS, ALPHA, TEMP, J, n_theta, 1009, 0.6, E, A, "full")
                r.update(pdb=pdb, sigma=g["sigma"], lam_max=g["lam_max"])
                rows.append(r)
            r = measure_hop(n, LAYERS, ALPHA, TEMP, 1.0, n_theta, 1009, 0.6, E, A, "hop_only")
            r.update(pdb=pdb, sigma=g["sigma"], lam_max=g["lam_max"]); rows.append(r)
            r = measure_hop(n, LAYERS, 1.0, 0.0, 0.0, n_theta, 1009, 0.6, E, A, "linear")
            r.update(pdb=pdb, sigma=g["sigma"], lam_max=g["lam_max"]); rows.append(r)
        print(f"  {pdb} done ({(time.time()-t0)/60:.1f} min)", flush=True)
    # summary: median over targets per (mode, J, n); fitted log2 slope per qubit
    summ = {}
    cells = sorted({(r["mode"], r["J"]) for r in rows})
    for mode, J in cells:
        per_n = {}
        for n in ns:
            v = [r["var_g0"] for r in rows if r["mode"] == mode and r["J"] == J and r["n"] == n]
            m = [r["mean_sq_per_param"] for r in rows if r["mode"] == mode and r["J"] == J and r["n"] == n]
            per_n[str(n)] = dict(var_g0_median=float(np.median(v)), var_g0_min=float(min(v)),
                                 var_g0_max=float(max(v)), msq_median=float(np.median(m)))
        ns_ = np.array(ns, float)
        v = np.array([per_n[str(n)]["var_g0_median"] for n in ns])
        ok = v > 0
        slope = float(np.polyfit(ns_[ok], np.log2(v[ok]), 1)[0]) if ok.sum() > 2 else float("nan")
        summ[f"{mode}|J{J:g}"] = dict(per_n=per_n, log2_slope_per_qubit=slope)
    base = summ["full|J0"]["per_n"]
    for kx, s in summ.items():
        s["ratio_to_J0_per_n"] = {n: s["per_n"][n]["var_g0_median"] / base[n]["var_g0_median"] for n in base}
    out = dict(kind="property measurement, gradient variance vs J, no RMSD", lane="S28B",
               targets=pick, n_theta=n_theta, alpha=ALPHA, T=TEMP, layers=LAYERS, ns=ns,
               j_grid=list(J_GRID), rows=rows, summary=summ)
    from s24 import stats_lib as ST
    ST.save_atomic(TRAIN_OUT, out, module_file=__file__)
    print(f"\n{'cell':16s} " + " ".join(f"{'n='+str(n):>11}" for n in ns) + f" {'slope':>8}")
    for kx, s in summ.items():
        print(f"{kx:16s} " + " ".join(f"{s['per_n'][str(n)]['var_g0_median']:11.3e}" for n in ns)
              + f" {s['log2_slope_per_qubit']:+8.3f}")
    print("wrote", TRAIN_OUT)


# ================================================================================ drivers
def _done_pdbs(path: str) -> set:
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["pdb"])
                except Exception:
                    pass
    return done


def run_main(limit: int = 0, pdbs_only: Optional[List[str]] = None):
    from s25 import phys_lib as P
    os.makedirs(RESULTS, exist_ok=True)
    done = _done_pdbs(ROWS)
    pdbs = P.targets()
    if pdbs_only:
        pdbs = [p for p in pdbs if p in set(pdbs_only)]
    if limit:
        pdbs = pdbs[:limit]
    t0 = time.time()
    for i, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        t1 = time.time()
        rows = run_target(pdb)
        with open(ROWS + ".tmp", "w", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        with open(ROWS, "a", encoding="utf-8") as fh:
            with open(ROWS + ".tmp", encoding="utf-8") as src:
                fh.write(src.read())
        os.remove(ROWS + ".tmp")
        j0 = [r for r in rows if r["source"] == "vqe" and r["J"] == 0.0 and r["seed"] == 0][0]
        print(f"  [{i+1}/{len(pdbs)}] {pdb} n={rows[0]['n']} rows={len(rows)} "
              f"vqe0={j0['R1']['rmsd']:.3f} dis75={j0['rmsd_dis75']:.3f} m0={j0['m']} "
              f"{time.time()-t1:.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", ROWS)


def load_rows(path: str = ROWS) -> List[Dict]:
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def chain_main(arms: List[str]):
    """Built chain for named arms `source|s<seed>|graph|J<j>|R<k>`; resumable per (arm, pdb)."""
    from s25 import phys_lib as P
    rows = load_rows()
    by = {(arm_of(r), r["pdb"]): r for r in rows}
    done = set()
    if os.path.exists(CHAIN_ROWS):
        for r in load_rows(CHAIN_ROWS):
            done.add((r["arm"], r["pdb"]))
    pdbs = P.targets()
    t0 = time.time()
    for i, pdb in enumerate(pdbs):
        todo = [a for a in arms if (a, pdb) not in done]
        if not todo:
            continue
        cand = H.Candidates.from_universe(pdb, k=500)
        out = []
        for a in todo:
            base, R = a.rsplit("|", 1)
            r = by[(base, pdb)]
            C = np.asarray(r[R]["C"], float)
            t1 = time.time()
            ca = H.readout_projected(cand, C)
            out.append(dict(arm=a, pdb=pdb, n=int(cand.n), fold=int(cand.fold), basis="built_chain",
                            label="ACHIEVABLE", rmsd_cloud=oracle_rmsd(cand, C),
                            rmsd_chain=oracle_rmsd(cand, ca), fail18=bool(pdb in FAIL18),
                            secs=float(time.time() - t1)))
        with open(CHAIN_ROWS, "a", encoding="utf-8") as fh:
            for r in out:
                fh.write(json.dumps(r) + "\n")
        print(f"  [chain {i+1}/{len(pdbs)}] {pdb} {len(out)} arms "
              f"{sum(r['secs'] for r in out):.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", CHAIN_ROWS)


# ================================================================================ analysis
def analyse_main():
    from s24 import stats_lib as ST
    rows = load_rows()
    by = {}
    for r in rows:
        by.setdefault(arm_of(r), {})[r["pdb"]] = r
    pdbs = sorted(by["vqe|s0|NONE|J0"])
    n = len(pdbs)
    folds = ST.pinned_folds(pdbs)
    fail = np.array([p in FAIL18 for p in pdbs])

    def vec(arm, R):
        return np.array([by[arm][p][R]["rmsd"] for p in pdbs])

    def compact(o):
        return dict(effect=o["effect"], se=o["se"], mde=o["mde"], x_mde=o["effect_over_mde"],
                    ci_iid=o["ci95_iid"], ci_fold=o["ci95_fold"], folds_same_sign=o["folds_same_sign"],
                    W=o["n_better"], L=o["n_worse"], T=o["n_tied"], verdict=o["verdict"],
                    mean_a=o["mean_a"], mean_b=o["mean_b"], median_effect=o["median_effect"],
                    fmt=ST.fmt(o))

    out = dict(n=n, arms={}, contrasts={}, anchors={})
    # anchors
    out["anchors"]["dis75_mean"] = float(np.mean([by["vqe|s0|NONE|J0"][p]["rmsd_dis75"] for p in pdbs]))
    for s in SEEDS:
        out["anchors"][f"vqe_J0_R1_s{s}"] = float(vec(f"vqe|s{s}|NONE|J0", "R1").mean())
    # per-arm means and native-free diagnostics
    for arm, d in by.items():
        if len(d) != n:
            continue
        rs = [d[p] for p in pdbs]
        out["arms"][arm] = dict(
            n=n, **{f"mean_{R}": float(vec(arm, R).mean()) for R in ("R1", "R2", "R3")},
            **{f"fail18_{R}": float(vec(arm, R)[fail].mean()) for R in ("R1", "R2", "R3")},
            **{f"other_{R}": float(vec(arm, R)[~fail].mean()) for R in ("R1", "R2", "R3")},
            m_mean=float(np.mean([r["m"] for r in rs])), pr_median=float(np.median([r["pr"] for r in rs])),
            jac75_mean=float(np.mean([r["jac75"] for r in rs])),
            mass_top75_mean=float(np.mean([r["mass_top75"] for r in rs])),
            entropy_bits_mean=float(np.mean([r["entropy_bits"] for r in rs])),
            gate_pass_all=bool(all(r["gate_pass"] for r in rs)),
            gate_equality_frac=float(np.mean([r["gate_equality"] for r in rs])),
            holes_max=int(max(r["gate_n_holes"] for r in rs)),
            tv_vs_J0_mean=(float(np.mean([r["tv_vs_J0"] for r in rs])) if rs[0].get("tv_vs_J0") is not None else None),
            hop_mean=(float(np.mean([r["hop"] for r in rs])) if "hop" in rs[0] else None),
            hop_abs_mean=(float(np.mean([r["hop_abs"] for r in rs])) if "hop_abs" in rs[0] else None),
            sign_coh_mean=(float(np.mean([r["sign_coh"] for r in rs])) if "sign_coh" in rs[0] else None),
            R1_mean_pair_rmsd=float(np.mean([r["R1"]["mean_pair_rmsd"] for r in rs])),
            R3_mean_pair_rmsd=float(np.mean([r["R3"]["mean_pair_rmsd"] for r in rs])),
            R1_overlap_dis75=float(np.mean([r["R1"]["overlap_dis75"] for r in rs])),
            R3_overlap_dis75=float(np.mean([r["R3"]["overlap_dis75"] for r in rs])),
            secs_mean=float(np.mean([r["secs"] for r in rs])),
            pad_mass_max=float(max(r["pad_mass"] for r in rs)),
        )
    # F1: each VQE hopping arm vs its J=0 counterpart, same readout and seed
    for s in SEEDS:
        for g in GRAPHS:
            for J in J_GRID[1:]:
                arm = f"vqe|s{s}|{g}|J{J:g}"
                if arm not in out["arms"]:
                    continue
                for R in ("R1", "R2", "R3"):
                    o = ST.compare(vec(arm, R), vec(f"vqe|s{s}|NONE|J0", R), folds=folds, names=pdbs,
                                   label=f"F1 {arm} {R} - vqe|s{s}|NONE|J0 {R} (point cloud)")
                    out["contrasts"][f"F1|{arm}|{R}"] = compact(o)
                    # the readout question: every arm against the DEPLOYED J=0 R1
                    if R != "R1":
                        o = ST.compare(vec(arm, R), vec(f"vqe|s{s}|NONE|J0", "R1"), folds=folds, names=pdbs,
                                       label=f"F1b {arm} {R} - vqe|s{s}|NONE|J0 R1 (point cloud)")
                        out["contrasts"][f"F1b|{arm}|{R}"] = compact(o)
                # F3 / F4: REAL vs PERM, REAL vs RAND at the same J
                if g == "REAL":
                    for ctrl in ("PERM", "RAND"):
                        carm = f"vqe|s{s}|{ctrl}|J{J:g}"
                        if carm in out["arms"]:
                            for R in ("R1", "R2", "R3"):
                                o = ST.compare(vec(arm, R), vec(carm, R), folds=folds, names=pdbs,
                                               label=f"F3/4 {arm} {R} - {carm} {R} (point cloud)")
                                out["contrasts"][f"CTRL|{arm}|{ctrl}|{R}"] = compact(o)
                    # F2: VQE vs GS at the same J, R2 and R3
                    garm = f"gs|s-1|{g}|J{J:g}"
                    if garm in out["arms"]:
                        for R in ("R2", "R3"):
                            o = ST.compare(vec(arm, R), vec(garm, R), folds=folds, names=pdbs,
                                           label=f"F2 {arm} {R} - {garm} {R} (point cloud)")
                            out["contrasts"][f"F2|{arm}|{R}"] = compact(o)
    # GS arms against the deployed J=0 VQE R1 (the eigensolver as a selector)
    for g in GRAPHS:
        for J in J_GRID[1:]:
            garm = f"gs|s-1|{g}|J{J:g}"
            if garm in out["arms"]:
                for R in ("R1", "R2", "R3"):
                    o = ST.compare(vec(garm, R), vec("vqe|s0|NONE|J0", "R1"), folds=folds, names=pdbs,
                                   label=f"GS {garm} {R} - vqe|s0|NONE|J0 R1 (point cloud)")
                    out["contrasts"][f"GS|{garm}|{R}"] = compact(o)
    # seed replication: seed-1 effect inside seed-0's fold CI
    for key, c in list(out["contrasts"].items()):
        if key.startswith("F1|vqe|s0|"):
            k1 = key.replace("|s0|", "|s1|")
            if k1 in out["contrasts"]:
                e1 = out["contrasts"][k1]["effect"]
                ci0 = c["ci_fold"]
                c["replicates_seed1"] = bool(ci0[0] <= e1 <= ci0[1])
                c["effect_seed1"] = e1
    # the J grid priced as an order statistic (REAL, R1 and R3, seed 0): min over J
    for R in ("R1", "R3"):
        M = np.stack([vec(f"vqe|s0|REAL|J{J:g}", R) - vec("vqe|s0|NONE|J0", R) for J in J_GRID[1:]], 1)
        out["contrasts"][f"GRID|REAL|s0|{R}"] = ST.best_of_k_within(M)
    # chain rows if present
    if os.path.exists(CHAIN_ROWS):
        cr = load_rows(CHAIN_ROWS)
        byc = {}
        for r in cr:
            byc.setdefault(r["arm"], {})[r["pdb"]] = r
        out["chain"] = {}
        pc = sorted(next(iter(byc.values())))
        pc = [p for p in pc if all(p in d for d in byc.values())]
        fc = ST.pinned_folds(pc)
        failc = np.array([p in FAIL18 for p in pc])
        for arm, d in byc.items():
            if not all(p in d for p in pc):
                continue
            a = np.array([d[p]["rmsd_chain"] for p in pc]); k = np.array([d[p]["rmsd_cloud"] for p in pc])
            out["chain"][arm] = dict(n=len(pc), mean_chain=float(a.mean()), mean_cloud=float(k.mean()),
                                     fail18_chain=float(a[failc].mean()), other_chain=float(a[~failc].mean()))
        ref = "vqe|s0|NONE|J0|R1"
        if ref in byc:
            base = np.array([byc[ref][p]["rmsd_chain"] for p in pc])
            for arm, d in byc.items():
                if arm == ref or not all(p in d for p in pc):
                    continue
                a = np.array([d[p]["rmsd_chain"] for p in pc])
                comp = ref.replace("|s0|", "|s1|") if "|s1|" in arm else ref
                b = base if comp == ref else np.array([byc[comp][p]["rmsd_chain"] for p in pc]) if comp in byc else base
                o = ST.compare(a, b, folds=fc, names=pc, label=f"CHAIN {arm} - {comp} (built chain)")
                out["chain"][arm]["vs_J0"] = compact(o)
                if failc.sum() > 2:
                    o2 = ST.compare(a[~failc], b[~failc], folds=fc[~failc], names=[p for p, f in zip(pc, failc) if not f],
                                    label=f"CHAIN {arm} - {comp} (built chain, 108 non-FAIL18)")
                    out["chain"][arm]["vs_J0_nonfail18"] = compact(o2)
    ST.save_atomic(SUMMARY, out, module_file=__file__)
    print(f"n={n}  DIS top-75 {out['anchors']['dis75_mean']:.6f}  VQE J0 R1 s0 {out['anchors']['vqe_J0_R1_s0']:.6f}  s1 {out['anchors']['vqe_J0_R1_s1']:.6f}")
    print(f"\n{'arm':22s} {'R1':>7} {'R2':>7} {'R3':>7} {'m':>5} {'PR':>7} {'jac75':>6} {'mass75':>6} {'Hbits':>6} {'TV':>6} {'hop':>6} {'coh':>5} {'eq':>4}")
    for arm in sorted(out["arms"], key=lambda a: (a.split("|")[0], a.split("|")[2], float(a.split("|")[3][1:]), a.split("|")[1])):
        d = out["arms"][arm]
        print(f"{arm:22s} {d['mean_R1']:7.4f} {d['mean_R2']:7.4f} {d['mean_R3']:7.4f} {d['m_mean']:5.1f} {d['pr_median']:7.1f} "
              f"{d['jac75_mean']:6.3f} {d['mass_top75_mean']:6.3f} {d['entropy_bits_mean']:6.2f} "
              f"{(d['tv_vs_J0_mean'] if d['tv_vs_J0_mean'] is not None else float('nan')):6.3f} "
              f"{(d['hop_mean'] if d['hop_mean'] is not None else float('nan')):6.3f} "
              f"{(d['sign_coh_mean'] if d['sign_coh_mean'] is not None else float('nan')):5.2f} {d['gate_equality_frac']:4.2f}")
    print("\nF1 (VQE hopping arm vs its own J=0, same readout and seed), point cloud:")
    for key, c in out["contrasts"].items():
        if key.startswith("F1|"):
            print(f"  {key:34s} eff {c['effect']:+.4f} x{c['x_mde']:+.2f} fold {c['ci_fold'][0]:+.4f},{c['ci_fold'][1]:+.4f} "
                  f"{c['W']}W/{c['L']}L/{c['T']}T rep={c.get('replicates_seed1')}  {c['verdict']}")
    print("wrote", SUMMARY)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--train", action="store_true")
    ap.add_argument("--chain", action="store_true")
    ap.add_argument("--analyse", action="store_true")
    ap.add_argument("--arms", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--pdbs", default="")
    ap.add_argument("--n-theta", type=int, default=120)
    a = ap.parse_args()
    if a.run:
        run_main(limit=a.limit, pdbs_only=[p for p in a.pdbs.split(",") if p] or None)
    if a.train:
        train_main(n_theta=a.n_theta)
    if a.chain:
        chain_main([x for x in a.arms.split(",") if x])
    if a.analyse:
        analyse_main()


if __name__ == "__main__":
    main()
