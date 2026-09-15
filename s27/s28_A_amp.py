#!/usr/bin/env python
"""s27/s28_A_amp.py -- S28 lane A: the AMPLITUDE READOUT.

A CVaR-VQE whose real-amplitude state is read as a SIGNED (affine) combination of the pool
windows,   C(theta) = sum_i psi_i(theta) W_i / sum_i psi_i(theta),   optimised under
    F_lam(theta) = CVaR_alpha(E; psi^2) - T H(psi^2) + lam * S~(C(theta)),
with S~ the shipped distogram Bayes risk of the emitted structure's own distance map (linearly
interpolated risk table, ADDENDUM 1 of `s27/PREREG_S28_A.md`).  Two halves, measured apart and
labelled apart: EXPRESSIVITY (ORACLE: the family's ceiling, theta chosen against the native) and
RECOGNITION (deployable: the native-free objective).  Classical controls replace ONLY the map
theta -> psi and share every other line of code.

Pre-registration: `s27/PREREG_S28_A.md` (read it first).  Tests: `tests/test_s28_A.py`.

Phases (each a governed job, per-target checkpoints, resumable):
    python s27/s28_A_amp.py oracle  [--limit N] [--reverse]   # frame, production, ORACLE ceilings
    python s27/s28_A_amp.py recog   [--limit N] [--seed S]    # every recognition arm, point cloud
    python s27/s28_A_amp.py chain   [--arms a,b,c]            # built chain of stored structures
    python s27/s28_A_amp.py analyse                           # ST.fmt blocks + summary json

ORACLE reads of `nat_ca` happen only inside functions whose names begin `oracle_`.  Everything
else is native-free, and the NaN-poison test in `tests/test_s28_A.py` asserts it.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q              # noqa: E402
from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s22 import qcand_lib as QC            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
STRUCTS = os.path.join(RESULTS, "s28_A_structs")
SALT = "s28A"
K = 500
M = 75
N_QUBITS, LAYERS = 9, 3
ALPHA, TEMP = 0.18, 0.5
LR, ITERS, ITERS_LONG = 0.15, 80, 400
LAMS = (0.3, 1.0, 3.0)
N_SUB, SUB_DIM = 8, 27
N_UNTRAINED = 16
ORACLE_STARTS, ORACLE_ITERS, ORACLE_LR = 5, 300, 0.05
BUDGET_MATCHED, BUDGET_CONV = 80, 2000
DENOM_EPS = 1e-9
#: the primary built-chain list (PREREG section 4, job 2, unconditional)
CHAIN_PRIMARY = (["prod", "circ_l0_i80"] + ["circ_l%g_i80" % l for l in LAMS]
                 + ["a500_rand_m_l%g" % l for l in LAMS] + ["a75_rand_m_l%g" % l for l in LAMS]
                 + ["untr_0", "oracle_circ", "oracle_aff500"])


# ======================================================================== pool and frame
def load_pool(pdb):
    """The shipped pool, its DIS score (asserted against the S25 cache) and the distogram."""
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(pdb)
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    top = RP.topm(ch["DIS"], M, key)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    return cand, np.asarray(ch["DIS"], float), np.sort(top).astype(int), dg


def zrank(x):
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(float(r.std()), 1e-12)


class Frame:
    """All k windows posed on the medoid of the retained set (the deployed average's frame).

    NATIVE-FREE: `top` is chosen by the score, the medoid by the set.  `Wp` (k, n, 3) are the
    posed windows, `Wf` (k, 3n) their flattening, `b` the medoid's index in `top`.
    """

    def __init__(self, W, top):
        W = np.asarray(W, float)
        top = np.asarray(top, int)
        P = I.pairwise_rmsd(W[top])
        self.b = int(I.medoid(P))
        self.ref = W[top][self.b]
        self.Wp = I.superpose_batch(W, self.ref)
        self.k, self.n = int(self.Wp.shape[0]), int(self.Wp.shape[1])
        self.Wf = self.Wp.reshape(self.k, 3 * self.n)
        self.top = top


def readout(psi, frame: Frame):
    """The amplitude readout: affine weights w = psi / sum psi over the REAL candidates.

    `psi` is (k,) real (padding already stripped).  Returns (C (n,3), denom, w).  If the
    denominator is below DENOM_EPS the readout is undefined: C is NaN and the caller records it.
    """
    psi = np.asarray(psi, float)
    denom = float(psi.sum())
    if abs(denom) < DENOM_EPS:
        return np.full((frame.n, 3), np.nan), denom, np.full(frame.k, np.nan)
    w = psi / denom
    C = (w @ frame.Wf).reshape(frame.n, 3)
    return C, denom, w


def weight_diag(w):
    w = np.asarray(w, float)
    if not np.isfinite(w).all():
        return dict(frac_neg=float("nan"), neg_mass=float("nan"), ess=float("nan"))
    return dict(frac_neg=float((w < 0).mean()), neg_mass=float(-w[w < 0].sum()),
                ess=float(1.0 / max((w ** 2).sum(), 1e-300)))


def struct_diag(C):
    C = np.asarray(C, float)
    if not np.isfinite(C).all():
        return dict(rg=float("nan"), bond=float("nan"))
    rg = float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean()))
    bond = float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean())
    return dict(rg=rg, bond=bond)


# ======================================================================= the objective term
class Surrogate:
    """S~(C): the shipped risk table read by linear interpolation, extended linearly beyond the
    grid ends (PREREG ADDENDUM 1).  `shipped(C)` is the exact shipped lookup, for reporting."""

    def __init__(self, dg, n):
        self.risk = np.asarray(dg["risk"], float)            # (npairs, G)
        self.grid = np.asarray(dg["grid"], float)
        self.g0, self.dg = float(self.grid[0]), float(self.grid[1] - self.grid[0])
        self.G = int(self.risk.shape[1])
        self.i, self.j = I.pair_index(int(n))
        assert np.array_equal(np.asarray(dg["i"], int), self.i) and \
            np.array_equal(np.asarray(dg["j"], int), self.j), "distogram pair order mismatch"
        self.npairs = int(len(self.i))
        self.slope = np.diff(self.risk, axis=1) / self.dg  # (npairs, G-1)
        self.dgd = dg

    def value_grad(self, C):
        """(S~, dS~/dC (n,3)).  Piecewise linear in every pair distance."""
        C = np.asarray(C, float)
        v = C[self.i] - C[self.j]
        d = np.linalg.norm(v, axis=1)
        u = (d - self.g0) / self.dg
        i0 = np.clip(np.floor(u).astype(int), 0, self.G - 2)
        f = u - i0                                         # may be <0 or >1 beyond the grid
        r0 = self.risk[np.arange(self.npairs), i0]
        s = self.slope[np.arange(self.npairs), i0]
        val = float((r0 + s * f * self.dg).mean())
        dS_dd = s / self.npairs
        dd_dv = v / np.maximum(d, 1e-12)[:, None]
        g = np.zeros_like(C)
        np.add.at(g, self.i, dS_dd[:, None] * dd_dv)
        np.add.at(g, self.j, -dS_dd[:, None] * dd_dv)
        return val, g

    def shipped(self, C):
        """The shipped score exactly as `s24.d_harness.score_shipped` computes it."""
        C = np.asarray(C, float)
        if not np.isfinite(C).all():
            return float("nan")
        D = I.pair_dists(C[None], self.i, self.j)
        return float(I.shipped_score(self.dgd, D.astype(np.float32).astype(float))[0])


# ============================================================ the objective, two families
def cvar_entropy_terms(E, p, alpha, T):
    """CVaR_alpha(E; p) - T H(p) and its gradient in p (the same formula as `core.quantum.free_energy`)."""
    v, q, dp = Q.cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-15))
    H = float(-(p * lp).sum())
    dH = -(lp + 1.0)
    return float(v - T * H), dp - T * dH, float(v), H


def objective_theta(circ, theta, E, alpha, T, lam, frame: Frame, sur: Surrogate):
    """F_lam(theta) and its EXACT gradient.  Returns (F, grad, info)."""
    f_fe, g_fe, p, cv, Hn = Q.free_energy(circ, theta, E, alpha, T)
    info = dict(cvar=float(cv), H=float(Hn), pad_mass=float(p[frame.k:].sum()))
    if lam == 0.0:
        return float(f_fe), g_fe, info
    psi = circ.state(theta)
    C, denom, w = readout(psi[:frame.k], frame)
    if not np.isfinite(C).all():
        info.update(S_smooth=float("nan"), denom=denom, undefined=True)
        return float(f_fe) + lam * 1e3, g_fe, info
    s_val, dS_dC = sur.value_grad(C)
    gC = dS_dC.ravel()
    dS_dpsi = (frame.Wf @ gC - float(C.ravel() @ gC)) / denom
    P = theta.size
    J = 0.5 * circ.states_batch(theta[None, :] + math.pi * np.eye(P))    # (P, dim) = dpsi/dtheta
    g_S = J[:, :frame.k] @ dS_dpsi
    info.update(S_smooth=float(s_val), denom=denom, undefined=False)
    return float(f_fe) + lam * s_val, g_fe + lam * g_S, info


def objective_psi(psi, E, alpha, T, lam, frame: Frame, sur: Surrogate):
    """The SAME objective for an arbitrary real amplitude vector psi (dim = len(E)).

    p = psi^2/|psi|^2; the readout uses psi[:k].  Returns (F, dF/dpsi, info).  This is the code
    path every classical control shares; a control differs from the circuit only in how psi is
    produced from its parameters.
    """
    psi = np.asarray(psi, float)
    n2 = float(psi @ psi)
    p = psi ** 2 / n2
    f_fe, gp, cv, Hn = cvar_entropy_terms(E, p, alpha, T)
    g = 2.0 * psi * (gp - float(gp @ p)) / n2
    info = dict(cvar=cv, H=Hn, pad_mass=float(p[frame.k:].sum()))
    if lam == 0.0:
        return f_fe, g, info
    C, denom, w = readout(psi[:frame.k], frame)
    if not np.isfinite(C).all():
        info.update(S_smooth=float("nan"), denom=denom, undefined=True)
        return f_fe + lam * 1e3, g, info
    s_val, dS_dC = sur.value_grad(C)
    gC = dS_dC.ravel()
    dS_dpsi = (frame.Wf @ gC - float(C.ravel() @ gC)) / denom
    g = g.copy()
    g[:frame.k] += lam * dS_dpsi
    info.update(S_smooth=float(s_val), denom=denom, undefined=False)
    return f_fe + lam * s_val, g, info


def objective_w(w, frame: Frame, sur: Surrogate):
    """The brief's literal control: S~ alone over affine weights w (k,), C = sum w W / sum w."""
    w = np.asarray(w, float)
    C, denom, ww = readout(w, frame)
    if not np.isfinite(C).all():
        return 1e3, np.zeros_like(w), dict(S_smooth=float("nan"), denom=denom, undefined=True)
    s_val, dS_dC = sur.value_grad(C)
    gC = dS_dC.ravel()
    g = (frame.Wf @ gC - float(C.ravel() @ gC)) / denom
    return s_val, g, dict(S_smooth=float(s_val), denom=denom, undefined=False)


# ============================================================================ optimisers
def adam(fg, x0, iters, lr, trace_every=0):
    """The deployed Adam loop of `core.quantum.run_cvar_vqe`, verbatim in its arithmetic."""
    x = np.asarray(x0, float).copy()
    m = np.zeros_like(x)
    v = np.zeros_like(x)
    trace = []
    for t in range(1, iters + 1):
        f, g, _ = fg(x)
        if trace_every and (t - 1) % trace_every == 0:
            trace.append((t - 1, float(f)))
        m = 0.9 * m + 0.1 * g
        v = 0.999 * v + 0.001 * g * g
        x = x - lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
    f, g, info = fg(x)
    if trace_every:
        trace.append((iters, float(f)))
    return x, float(f), info, trace


def lbfgs(fg, x0, maxfun):
    """L-BFGS-B with a hard evaluation budget.  Returns (x, f, info, n_evals)."""
    from scipy.optimize import minimize
    cnt = [0]

    def fun(x):
        cnt[0] += 1
        f, g, _ = fg(x)
        return float(f), np.asarray(g, float)
    r = minimize(fun, np.asarray(x0, float), jac=True, method="L-BFGS-B",
                 options=dict(maxfun=int(maxfun), maxiter=int(maxfun), ftol=1e-12, gtol=1e-9))
    f, g, info = fg(r.x)
    return np.asarray(r.x, float), float(f), info, int(cnt[0])


# =============================================================== the classical families
class Family:
    """A map z -> psi (dim) with its chain rule, for the matched-objective controls.

    Every family lives on the k REAL candidates: psi[k:] (the padding amplitudes) is exactly
    zero, so the padding carries no probability and the objective is the circuit's, verbatim.
    """

    def __init__(self, kind, dim, k, top, B=None):
        self.kind, self.dim, self.k, self.top, self.B = kind, int(dim), int(k), top, B

    def n_params(self):
        if self.kind == "sub":
            return int(self.B.shape[1])
        return {"a500": self.k, "a75": len(self.top), "simplex": self.k}[self.kind]

    def psi(self, z):
        psi = np.zeros(self.dim)
        if self.kind == "a500":
            psi[:self.k] = z
        elif self.kind == "a75":
            psi[self.top] = z
        elif self.kind == "sub":
            psi[:self.k] = self.B @ np.asarray(z, float)
        elif self.kind == "simplex":
            z = np.asarray(z, float)
            e = np.exp(z - z.max())
            psi[:self.k] = np.sqrt(e / e.sum())
        else:
            raise KeyError(self.kind)
        return psi

    def chain(self, z, psi, g_psi):
        if self.kind == "a500":
            return g_psi[:self.k]
        if self.kind == "a75":
            return g_psi[self.top]
        if self.kind == "sub":
            return self.B.T @ g_psi[:self.k]
        if self.kind == "simplex":
            ps, gs = psi[:self.k], g_psi[:self.k]
            return 0.5 * (gs * ps - (ps ** 2) * float(gs @ ps))
        raise KeyError(self.kind)

    def init(self, which, rng):
        n = self.n_params()
        if which == "rand":
            return rng.normal(0.0, 1.0, n)
        if which == "prod":
            if self.kind == "a500":
                z = np.zeros(self.k); z[self.top] = 1.0; return z
            if self.kind == "a75":
                return np.ones(len(self.top))
            if self.kind == "simplex":
                z = np.full(self.k, -40.0); z[self.top] = 0.0; return z
            raise ValueError("a random subspace has no production point")
        raise KeyError(which)


class WFamily:
    """A map z -> affine weights w (k,) for the S-only controls (the brief's literal form)."""

    def __init__(self, kind, k, top, B=None):
        self.kind, self.k, self.top, self.B = kind, int(k), top, B

    def n_params(self):
        if self.kind == "sub":
            return int(self.B.shape[1])
        return {"a500": self.k, "a75": len(self.top), "simplex": self.k}[self.kind]

    def w(self, z):
        if self.kind == "a500":
            return np.asarray(z, float)
        if self.kind == "a75":
            w = np.zeros(self.k); w[self.top] = z; return w
        if self.kind == "sub":
            return self.B @ np.asarray(z, float)
        if self.kind == "simplex":
            z = np.asarray(z, float)
            e = np.exp(z - z.max())
            return e / e.sum()
        raise KeyError(self.kind)

    def chain(self, z, w, g_w):
        if self.kind == "a500":
            return g_w
        if self.kind == "a75":
            return g_w[self.top]
        if self.kind == "sub":
            return self.B.T @ g_w
        if self.kind == "simplex":
            return w * (g_w - float(g_w @ w))
        raise KeyError(self.kind)

    def init(self, which, rng):
        n = self.n_params()
        if which == "rand":
            return rng.normal(0.0, 1.0, n)
        if which == "prod":
            if self.kind == "a500":
                z = np.zeros(self.k); z[self.top] = 1.0 / len(self.top); return z
            if self.kind == "a75":
                return np.full(len(self.top), 1.0 / len(self.top))
            if self.kind == "simplex":
                z = np.full(self.k, -40.0); z[self.top] = 0.0; return z
            raise ValueError("a random subspace has no production point")
        raise KeyError(which)


def subspace_matrix(pdb, s, rows, cols=SUB_DIM):
    return SD.stable_rng(pdb, "s28A_subspace", int(s), int(rows), salt=SALT).normal(0.0, 1.0, (rows, cols))


# ================================================================================ ORACLE
def oracle_kabsch(C, nat):
    """ORACLE.  (rmsd, R, residual) with residual = R Cc - Nc; d rmsd / dC = residual @ R / (n rmsd)."""
    C = np.asarray(C, float); nat = np.asarray(nat, float)
    Cc = C - C.mean(0); Nc = nat - nat.mean(0)
    Hm = Cc.T @ Nc
    U, S, Vt = np.linalg.svd(Hm)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    res = Cc @ R.T - Nc
    rmsd = float(np.sqrt((res ** 2).sum() / len(C)))
    return rmsd, R, res


def oracle_rmsd_grad(C, nat):
    """ORACLE.  CA-RMSD after optimal superposition and its gradient in C (envelope theorem)."""
    rmsd, R, res = oracle_kabsch(C, nat)
    g = (res @ R) / (len(C) * max(rmsd, 1e-12))
    return rmsd, g


def oracle_circuit_ceiling(circ, frame: Frame, nat, starts=ORACLE_STARTS, iters=ORACLE_ITERS,
                           lr=ORACLE_LR):
    """ORACLE.  Minimise the point-cloud RMSD of C(theta) to the native over theta; best of `starts`."""
    P = circ.n_params()
    per = []
    best = None
    for s in range(int(starts)):
        th0 = np.random.default_rng(int(s)).normal(0.0, 0.6, P)

        def fg(th):
            psi = circ.state(th)
            C, denom, w = readout(psi[:frame.k], frame)
            if not np.isfinite(C).all():
                return 1e3, np.zeros(P), {}
            r, gC = oracle_rmsd_grad(C, nat)
            gC = gC.ravel()
            d_dpsi = (frame.Wf @ gC - float(C.ravel() @ gC)) / denom
            J = 0.5 * circ.states_batch(th[None, :] + math.pi * np.eye(P))
            return r, J[:, :frame.k] @ d_dpsi, {}
        th, f, _, _ = adam(fg, th0, iters, lr)
        psi = circ.state(th)
        C, denom, w = readout(psi[:frame.k], frame)
        r = float(I.ca_rmsd(C, nat)) if np.isfinite(C).all() else float("nan")
        per.append(r)
        if best is None or (np.isfinite(r) and r < best[0]):
            best = (r, C, th, denom, w)
    return dict(rmsd_cloud=best[0], per_start=per, rmsd_mean_starts=float(np.nanmean(per)),
                C=best[1], theta=best[2], denom=best[3], **weight_diag(best[4]), **struct_diag(best[1]))


def oracle_affine_ls(frame: Frame, nat, idx=None):
    """ORACLE.  Least squares over the AFFINE hull of the windows `idx` (all if None), in the frame,
    against the native posed on the frame's reference.  Returns (rmsd after re-superposition, C, w)."""
    idx = np.arange(frame.k) if idx is None else np.asarray(idx, int)
    Wf = frame.Wf[idx]
    natp = I.superpose_batch(np.asarray(nat, float)[None], frame.ref)[0].ravel()
    base = Wf[-1]
    A = (Wf[:-1] - base).T                              # (3n, m-1)
    y, *_ = np.linalg.lstsq(A, natp - base, rcond=None)
    w = np.empty(len(idx)); w[:-1] = y; w[-1] = 1.0 - y.sum()
    C = (w @ Wf).reshape(frame.n, 3)
    return float(I.ca_rmsd(C, nat)), C, w


def oracle_subspace_ls(frame: Frame, nat, B):
    """ORACLE.  Least squares over w = B z with sum w = 1 (a random 27-dim affine subspace)."""
    natp = I.superpose_batch(np.asarray(nat, float)[None], frame.ref)[0].ravel()
    c = B.sum(0)                                        # sum w = c . z = 1
    z0 = c / float(c @ c)
    Qm, _ = np.linalg.qr(np.column_stack([c, np.eye(len(c))]))
    N = Qm[:, 1:len(c)]                                  # null space of c^T
    A = frame.Wf.T @ B                                  # (3n, 27)
    y, *_ = np.linalg.lstsq(A @ N, natp - A @ z0, rcond=None)
    z = z0 + N @ y
    w = B @ z
    C = (w @ frame.Wf).reshape(frame.n, 3)
    return float(I.ca_rmsd(C, nat)), C, w


# ============================================================================ per target
def _row_for(name, C, sur: Surrogate, extra=None):
    r = dict(S_shipped=sur.shipped(C), **struct_diag(C))
    if extra:
        r.update(extra)
    return r


def run_oracle_target(pdb, seed=0):
    """Phase 1: frame, production, the lam = 0 soundness gate, ORACLE ceilings.  Stores structures."""
    t0 = time.time()
    cand, dis, top, dg = load_pool(pdb)
    frame = Frame(cand.W, top)
    sur = Surrogate(dg, cand.n)
    enc = QC.Encoding(zrank(dis))
    circ = Q.StatevectorCircuit(N_QUBITS, LAYERS)
    assert enc.dim == circ.dim and enc.K == frame.k
    structs = {}
    row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), k=int(cand.k), medoid=int(frame.b),
               seed=int(seed), arms={})
    # production: uniform on the top-75, and the deployed operator side by side
    psi_u = np.zeros(frame.k); psi_u[top] = 1.0
    Cp, _, wp = readout(psi_u, frame)
    from s24 import d_harness as H
    Cd, bd = H.readout_uniform(cand, top)
    row["prod_frame_max_dev"] = float(np.abs(Cp - Cd).max())
    row["arms"]["prod"] = _row_for("prod", Cp, sur, dict(S_smooth=sur.value_grad(Cp)[0],
                                                          rmsd_cloud=oracle_rmsd_of(Cp, cand),
                                                          rmsd_cloud_deployed=oracle_rmsd_of(Cd, cand)))
    structs["prod"] = Cp
    # pool descriptors (native-free) and ORACLE pool labels
    row["pool_rg_mean"] = float(np.mean([struct_diag(w)["rg"] for w in frame.Wp]))
    row["pool_bond_mean"] = float(np.mean([struct_diag(w)["bond"] for w in frame.Wp]))
    row["oracle_pool_best"] = float(np.min(cand.oracle_rr))
    row["oracle_top75_best"] = float(np.min(cand.oracle_rr[top]))
    # ORACLE ceilings
    t1 = time.time()
    oc = oracle_circuit_ceiling(circ, frame, cand.nat_ca)
    structs["oracle_circ"] = oc.pop("C"); oc.pop("theta")
    oc["secs"] = time.time() - t1
    oc["S_shipped"] = sur.shipped(structs["oracle_circ"])
    oc["S_smooth"] = float(sur.value_grad(structs["oracle_circ"])[0]) if np.isfinite(structs["oracle_circ"]).all() else float("nan")
    row["arms"]["oracle_circ"] = oc
    r5, C5, w5 = oracle_affine_ls(frame, cand.nat_ca)
    structs["oracle_aff500"] = C5
    row["arms"]["oracle_aff500"] = _row_for("oracle_aff500", C5, sur, dict(rmsd_cloud=r5, S_smooth=sur.value_grad(C5)[0], **weight_diag(w5)))
    r75, C75, w75 = oracle_affine_ls(frame, cand.nat_ca, top)
    structs["oracle_aff75"] = C75
    row["arms"]["oracle_aff75"] = _row_for("oracle_aff75", C75, sur, dict(rmsd_cloud=r75, S_smooth=sur.value_grad(C75)[0], **weight_diag(w75)))
    subs = []
    for s in range(N_SUB):
        B = subspace_matrix(pdb, s, frame.k)
        rs, Cs, ws = oracle_subspace_ls(frame, cand.nat_ca, B)
        subs.append(rs)
    row["arms"]["oracle_sub"] = dict(rmsd_cloud=float(np.mean(subs)), per_sub=subs)
    # native-free: the ORACLE structure's objective value is a diagnostic of the objective (§6)
    row["secs"] = time.time() - t0
    _save_structs(pdb, structs, cand)
    return row


def run_recog_target(pdb, seed=0, long_diag=True):
    """Phase 2: every recognition arm on the point cloud; stores the structures for the chain."""
    t0 = time.time()
    cand, dis, top, dg = load_pool(pdb)
    frame = Frame(cand.W, top)
    sur = Surrogate(dg, cand.n)
    enc = QC.Encoding(zrank(dis))
    E = enc.E
    circ = Q.StatevectorCircuit(N_QUBITS, LAYERS)
    P = circ.n_params()
    structs = {}
    row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), k=int(cand.k), seed=int(seed), arms={})

    def finish(name, C, denom, w, info, extra=None):
        r = dict(rmsd_cloud=oracle_rmsd_of(C, cand), S_shipped=sur.shipped(C), denom=float(denom),
                 undefined=bool(not np.isfinite(C).all()), **weight_diag(w), **struct_diag(C))
        r.update({k: v for k, v in info.items() if k in ("cvar", "H", "pad_mass", "S_smooth")})
        if extra:
            r.update(extra)
        row["arms"][name] = r
        structs[name] = C

    # --- (2), (3): the circuit
    th0 = np.random.default_rng(int(seed)).normal(0.0, 0.6, P)
    for lam in (0.0,) + tuple(LAMS):
        t1 = time.time()
        th, f, info, trace = adam(lambda x: objective_theta(circ, x, E, ALPHA, TEMP, lam, frame, sur),
                                  th0, ITERS, LR, trace_every=20)
        psi = circ.state(th)
        C, denom, w = readout(psi[:frame.k], frame)
        extra = dict(F=float(f), trace=trace, secs=time.time() - t1, iters=ITERS)
        if lam == 0.0:
            p_ref, cv_ref, H_ref, _ = Q.run_cvar_vqe(E, ALPHA, TEMP, n=N_QUBITS, layers=LAYERS,
                                                     iters=ITERS, restarts=1, seed=int(seed), lr=LR)
            extra["gate_lam0_max_abs_dp"] = float(np.abs(circ.probs(th) - p_ref).max())
            info["S_smooth"] = float(sur.value_grad(C)[0]) if np.isfinite(C).all() else float("nan")
        finish("circ_l%g_i%d" % (lam, ITERS), C, denom, w, info, extra)
    if long_diag:
        t1 = time.time()
        th, f, info, trace = adam(lambda x: objective_theta(circ, x, E, ALPHA, TEMP, 1.0, frame, sur),
                                  th0, ITERS_LONG, LR, trace_every=20)
        psi = circ.state(th)
        C, denom, w = readout(psi[:frame.k], frame)
        finish("diag_circ_l1_i%d" % ITERS_LONG, C, denom, w, info,
               dict(F=float(f), trace=trace, secs=time.time() - t1, iters=ITERS_LONG))

    # --- (4): classical controls, matched objective
    fams = [("a500", Family("a500", enc.dim, frame.k, top)), ("a75", Family("a75", enc.dim, frame.k, top)),
            ("simplex", Family("simplex", enc.dim, frame.k, top))]
    fams += [("sub%d" % s, Family("sub", enc.dim, frame.k, top, B=subspace_matrix(pdb, s, frame.k)))
             for s in range(N_SUB)]
    for fname, fam in fams:
        for init in ("rand", "prod"):
            if init == "prod" and fam.kind == "sub":
                continue
            z0 = fam.init(init, SD.stable_rng(pdb, "s28A_init", fname, init, int(seed), salt=SALT))
            for lam in LAMS:
                def fg(z, lam=lam):
                    psi = fam.psi(z)
                    f, g, info = objective_psi(psi, E, ALPHA, TEMP, lam, frame, sur)
                    return f, fam.chain(z, psi, g), info
                for bud, tag in ((BUDGET_MATCHED, "m"), (BUDGET_CONV, "c")):
                    t1 = time.time()
                    z, f, info, nev = lbfgs(fg, z0, bud)
                    psi = fam.psi(z)
                    C, denom, w = readout(psi[:frame.k], frame)
                    finish("%s_%s_%s_l%g" % (fname, init, tag, lam), C, denom, w, info,
                           dict(F=float(f), n_evals=nev, secs=time.time() - t1))
    # --- (4, S-only): the brief's literal control
    wfams = [("a500", WFamily("a500", frame.k, top)), ("a75", WFamily("a75", frame.k, top)),
             ("simplex", WFamily("simplex", frame.k, top))]
    wfams += [("sub%d" % s, WFamily("sub", frame.k, top, B=subspace_matrix(pdb, s, frame.k)))
              for s in range(N_SUB)]
    for fname, fam in wfams:
        for init in ("rand", "prod"):
            if init == "prod" and fam.kind == "sub":
                continue
            z0 = fam.init(init, SD.stable_rng(pdb, "s28A_winit", fname, init, int(seed), salt=SALT))

            def fg(z):
                w = fam.w(z)
                f, g, info = objective_w(w, frame, sur)
                return f, fam.chain(z, w, g), info
            for bud, tag in ((BUDGET_MATCHED, "m"), (BUDGET_CONV, "c")):
                t1 = time.time()
                z, f, info, nev = lbfgs(fg, z0, bud)
                w = fam.w(z)
                C, denom, ww = readout(w, frame)
                finish("sonly_%s_%s_%s" % (fname, init, tag), C, denom, ww, info,
                       dict(F=float(f), n_evals=nev, secs=time.time() - t1))
    # --- (5): the untrained circuit
    for d in range(N_UNTRAINED):
        th = SD.stable_rng(pdb, "s28A_untrained", int(d), int(seed), salt=SALT).normal(0.0, 0.6, P)
        psi = circ.state(th)
        C, denom, w = readout(psi[:frame.k], frame)
        p = psi ** 2
        f_fe, _, cv, Hn = cvar_entropy_terms(E, p, ALPHA, TEMP)
        info = dict(cvar=cv, H=Hn, pad_mass=float(p[frame.k:].sum()),
                    S_smooth=float(sur.value_grad(C)[0]) if np.isfinite(C).all() else float("nan"))
        finish("untr_%d" % d, C, denom, w, info)
    row["secs"] = time.time() - t0
    _save_structs(pdb, structs, cand, merge=True)
    return row


def oracle_rmsd_of(C, cand):
    """ORACLE.  Point-cloud CA-RMSD of an emitted structure; NaN if the structure is undefined."""
    if cand.nat_ca is None or not np.isfinite(np.asarray(C, float)).all()             or not np.isfinite(np.asarray(cand.nat_ca, float)).all():
        return float("nan")
    return float(I.ca_rmsd(C, cand.nat_ca))


def _save_structs(pdb, structs, cand, merge=True):
    os.makedirs(STRUCTS, exist_ok=True)
    f = os.path.join(STRUCTS, f"{pdb}.npz")
    old = {}
    if merge and os.path.exists(f):
        z = np.load(f)
        old = {k: z[k] for k in z.files}
    old.update({k: np.asarray(v, float) for k, v in structs.items()})
    old["seq"] = np.array(cand.seq); old["fold"] = np.array(int(cand.fold)); old["n"] = np.array(int(cand.n))
    tmp = f + f".tmp{os.getpid()}.npz"
    np.savez_compressed(tmp, **old)
    os.replace(tmp, f)


def run_chain_target(pdb, arms):
    """Phase 3: the BUILT CHAIN of stored structures through the production projection.  ORACLE
    scoring of the emitted chain against the native, post hoc."""
    z = np.load(os.path.join(STRUCTS, f"{pdb}.npz"))
    seq, fold = str(z["seq"]), int(z["fold"])
    u = I.load_univ(pdb)
    nat = u["nat_ca"]
    out = {}
    for a in arms:
        if a not in z.files:
            out[a] = dict(rmsd_chain=float("nan"), missing=True)
            continue
        C = np.asarray(z[a], float)
        if not np.isfinite(C).all():
            out[a] = dict(rmsd_chain=float("nan"), undefined=True)
            continue
        t1 = time.time()
        pr = I.project(C, seq, fold)
        ca = np.asarray(pr["ca"], float)
        out[a] = dict(rmsd_chain=float(I.ca_rmsd(ca, nat)), rmsd_cloud=float(I.ca_rmsd(C, nat)),
                      secs=time.time() - t1, **struct_diag(ca))
    return dict(pdb=pdb, fold=fold, arms=out)


# ================================================================================= driver
def _jsonl_done(path, key="pdb"):
    done = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line); done[r[key]] = r
                except Exception:
                    pass
    return done


def _append(path, row):
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)) + "\n")


def phase(mode, pdbs, seed=0, arms=None, reverse=False):
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"s28_A_{mode}_rows.jsonl" if seed == 0 else f"s28_A_{mode}_rows_seed{seed}.jsonl")
    done = _jsonl_done(path)
    if reverse:
        pdbs = list(reversed(pdbs))
    t0 = time.time()
    for n, pdb in enumerate(pdbs):
        if pdb in done and mode != "chain":
            continue
        if mode == "chain" and pdb in done and all(a in done[pdb]["arms"] for a in arms):
            continue
        t1 = time.time()
        if mode == "oracle":
            row = run_oracle_target(pdb, seed=seed)
            msg = "prod %.3f  ORACLE circ %.3f (starts %s)  aff500 %.3f aff75 %.3f sub %.3f  framedev %.1e" % (
                row["arms"]["prod"]["rmsd_cloud"], row["arms"]["oracle_circ"]["rmsd_cloud"],
                " ".join("%.2f" % v for v in row["arms"]["oracle_circ"]["per_start"]),
                row["arms"]["oracle_aff500"]["rmsd_cloud"], row["arms"]["oracle_aff75"]["rmsd_cloud"],
                row["arms"]["oracle_sub"]["rmsd_cloud"], row["prod_frame_max_dev"])
        elif mode == "recog":
            row = run_recog_target(pdb, seed=seed)
            a = row["arms"]
            msg = "circ l0 %.2f l0.3 %.2f l1 %.2f l3 %.2f | a500m %.2f a75m %.2f sub0m %.2f simplexm %.2f | untr0 %.2f | gate %.1e" % (
                a["circ_l0_i80"]["rmsd_cloud"], a["circ_l0.3_i80"]["rmsd_cloud"], a["circ_l1_i80"]["rmsd_cloud"],
                a["circ_l3_i80"]["rmsd_cloud"], a["a500_rand_m_l1"]["rmsd_cloud"], a["a75_rand_m_l1"]["rmsd_cloud"],
                a["sub0_rand_m_l1"]["rmsd_cloud"], a["simplex_rand_m_l1"]["rmsd_cloud"], a["untr_0"]["rmsd_cloud"],
                a["circ_l0_i80"]["gate_lam0_max_abs_dp"])
        else:
            todo = [a for a in arms if not (pdb in done and a in done[pdb]["arms"])]
            row = run_chain_target(pdb, todo)
            if pdb in done:
                merged = dict(done[pdb]); merged["arms"] = dict(done[pdb]["arms"]); merged["arms"].update(row["arms"])
                row = merged
            msg = " ".join("%s %.3f" % (a, row["arms"][a].get("rmsd_chain", float("nan"))) for a in todo)
        _append(path, row)
        done[pdb] = row
        print(f"  [{mode} {n+1}/{len(pdbs)}] {pdb} n={row.get('n', '?')} {msg} {time.time()-t1:.1f}s "
              f"(elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["oracle", "recog", "chain", "analyse"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--arms", default=",".join(CHAIN_PRIMARY))
    ap.add_argument("--reverse", action="store_true")
    ap.add_argument("--pdbs", default="")
    a = ap.parse_args()
    from s25 import phys_lib as P
    pdbs = P.targets()
    if a.pdbs:
        pdbs = [p for p in a.pdbs.split(",") if p]
    if a.limit:
        pdbs = pdbs[:a.limit]
    if a.mode == "analyse":
        from s27 import s28_A_analyse as AN
        AN.main()
        return
    phase(a.mode, pdbs, seed=a.seed, arms=[x for x in a.arms.split(",") if x], reverse=a.reverse)


if __name__ == "__main__":
    main()
