#!/usr/bin/env python
"""s27/s28_A2_local.py -- S28 lane A2: is the objective LOCALLY informative at the production point?

`s27/PREREG_S28_A.md` ADDENDUM 4; brief `s27/briefs/S28A2.md`.  At the production point cloud C0
(the DIS top-75 uniform average in its medoid frame):

  A2.1  ORACLE DIAGNOSTIC: the cosine between the steepest-descent direction of the shipped
        Bayes risk (-dS~/dC, analytic, rigid-body components removed) and the ORACLE direction
        to the native (u = native aligned onto C0, minus C0; rigid-body removed).  The same for
        the S27 CA channels (RG_LAW, EXVOL analytic; DISTPOT, CONTACT, ENV, CAGEO are step
        functions and get a SMOOTHED central difference at h = 0.5 A, labelled so).
  A2.2  DEPLOYABLE (native-free): C(e) = C0 - e g/rms(g) for e in {0.1, 0.3, 1.0} A; a random
        direction of the same displacement (8 draws); the same step restricted to the circuit
        family (from theta_P, the family's nearest point to C0).  Point cloud here; the built
        chain through `s28_A_amp.py chain` (which reads this module's `<pdb>_a2.npz`).

    python s27/s28_A2_local.py cosine [--limit N]      # A2.1, ORACLE diagnostic, 126 rows
    python s27/s28_A2_local.py ladder [--limit N]      # A2.2, point cloud + stored clouds

`nat_ca` is read only inside functions named `oracle_*`.  The deployable ladder is NaN-poison
tested in `tests/test_s28_A2.py`.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

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
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A_amp as A             # noqa: E402

RESULTS = A.RESULTS
SALT = "s28A2"
STEPS = (0.1, 0.3, 1.0)
N_RAND = 8
FD_H = 0.5
CHANNELS_ANALYTIC = ("DIS", "RG_LAW", "EXVOL")
CHANNELS_FD = ("DISTPOT", "CONTACT", "ENV", "CAGEO")


# ============================================================ rigid-body removal, cosine
def rigid_basis(C0):
    """Orthonormal basis (3n, 6) of translations and infinitesimal rotations about the centroid."""
    C0 = np.asarray(C0, float)
    n = len(C0)
    c = C0 - C0.mean(0)
    cols = []
    for ax in range(3):
        t = np.zeros((n, 3)); t[:, ax] = 1.0
        cols.append(t.ravel())
    for ax in range(3):
        w = np.zeros(3); w[ax] = 1.0
        cols.append(np.cross(w[None, :], c).ravel())
    B = np.column_stack(cols)
    Qm, _ = np.linalg.qr(B)
    return Qm


def remove_rigid(v, C0):
    """Project the rigid-body components out of a coordinate field v (n,3) at C0."""
    B = rigid_basis(C0)
    x = np.asarray(v, float).ravel()
    return (x - B @ (B.T @ x)).reshape(np.asarray(v).shape)


def rms(v):
    v = np.asarray(v, float)
    return float(np.sqrt((v ** 2).sum(-1).mean()))


def cosine(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-15 or nb < 1e-15:
        return float("nan")
    return float(a @ b / (na * nb))


# ============================================================ gradients of the channels
def grad_rg_law(C, n):
    """RG_LAW = |Rg - r0| / r0, r0 = 2.2 n^0.38; analytic gradient."""
    C = np.asarray(C, float)
    r0 = 2.2 * n ** 0.38
    c = C - C.mean(0)
    rg = float(np.sqrt((c ** 2).sum(1).mean()))
    drg = c / (len(C) * max(rg, 1e-12))
    return abs(rg - r0) / r0, np.sign(rg - r0) * drg / r0


def grad_exvol(C, cut=4.0):
    """EXVOL = sum_{|i-j|>=3} max(0, cut - d)^2; analytic gradient."""
    C = np.asarray(C, float)
    i, j = I.pair_index(len(C), 3)
    v = C[i] - C[j]
    d = np.linalg.norm(v, axis=1)
    pen = np.maximum(0.0, cut - d)
    val = float((pen ** 2).sum())
    coef = -2.0 * pen / np.maximum(d, 1e-12)
    g = np.zeros_like(C)
    np.add.at(g, i, coef[:, None] * v)
    np.add.at(g, j, -coef[:, None] * v)
    return val, g


class _FakeCand:
    def __init__(self, cand, W):
        self.seq, self.n, self.fold = cand.seq, cand.n, cand.fold
        self.W = np.asarray(W, float)
        self.k = len(self.W)
        self.PHI = np.zeros((self.k, self.n)); self.PSI = np.zeros((self.k, self.n))


def fd_channel_grads(cand, C0, h=FD_H):
    """Smoothed central differences of the S27 step-function channels at C0, one batched call
    per channel (the universe fit is done once per call).  Labelled 'smoothed FD'."""
    from s27 import ham_lib as HL
    from s27 import run_pool as RP
    C0 = np.asarray(C0, float)
    n = len(C0)
    stack = [C0]
    for k in range(3 * n):
        e = np.zeros(3 * n); e[k] = h
        stack.append(C0 + e.reshape(n, 3)); stack.append(C0 - e.reshape(n, 3))
    W = np.stack(stack)
    u = I.load_univ(cand.pdb)
    dg = I.distogram(cand.pdb, cand.seq, cand.fold)
    cx = HL.Context(_FakeCand(cand, W), u, dg, RP.RAMA[cand.fold])
    fns = dict(DISTPOT=HL.h_distpot, CONTACT=HL.h_contact, ENV=HL.h_env, CAGEO=HL.h_cageo)
    out = {}
    for name in CHANNELS_FD:
        vals = np.asarray(fns[name](cx), float)
        g = (vals[1::2] - vals[2::2]) / (2 * h)
        out[name] = (float(vals[0]), g.reshape(n, 3))
    return out


# ============================================================ ORACLE pieces
def oracle_direction(C0, nat):
    """ORACLE.  u = (native Kabsch-aligned onto C0) - C0."""
    natp = I.superpose_batch(np.asarray(nat, float)[None], np.asarray(C0, float))[0]
    return natp - np.asarray(C0, float)


def oracle_rmsd(C, cand):
    return A.oracle_rmsd_of(C, cand)


# ============================================================ A2.1 per target
def production_point(pdb):
    cand, dis, top, dg = A.load_pool(pdb)
    frame = A.Frame(cand.W, top)
    sur = A.Surrogate(dg, cand.n)
    psi = np.zeros(frame.k); psi[top] = 1.0
    C0, _, _ = A.readout(psi, frame)
    return cand, frame, sur, C0, top


def run_cosine_target(pdb, with_fd=True):
    t0 = time.time()
    cand, frame, sur, C0, top = production_point(pdb)
    grads = {}
    s_val, g = sur.value_grad(C0)
    grads["DIS"] = (s_val, g)
    grads["RG_LAW"] = grad_rg_law(C0, cand.n)
    grads["EXVOL"] = grad_exvol(C0)
    if with_fd:
        grads.update(fd_channel_grads(cand, C0))
    u = remove_rigid(oracle_direction(C0, cand.nat_ca), C0)                 # ORACLE
    row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), fail18=bool(pdb in I.FAIL18),
               rmsd_prod_cloud=oracle_rmsd(C0, cand), u_rms=rms(u), cos={}, grad_rms={}, value={})
    for name, (val, gg) in grads.items():
        gp = remove_rigid(gg, C0)
        row["cos"][name] = cosine(-gp, u)                                    # ORACLE (u)
        row["grad_rms"][name] = rms(gp)
        row["value"][name] = float(val)
    # a random-direction reference for the cosine distribution (what |cos| a random field gives)
    rng = SD.stable_rng(pdb, "s28A2_cosref", salt=SALT)
    row["cos_random_ref"] = [cosine(remove_rigid(rng.normal(size=C0.shape), C0), u) for _ in range(16)]
    row["secs"] = time.time() - t0
    return row


# ============================================================ A2.2 per target (deployable)
def theta_nearest_production(circ, frame, C0, starts=5, iters=300, lr=0.05):
    """NATIVE-FREE: the circuit family's nearest point to the production cloud C0."""
    P = circ.n_params()
    best = None
    for s in range(starts):
        th0 = np.random.default_rng(int(s)).normal(0.0, 0.6, P)

        def fg(th):
            psi = circ.state(th)
            C, denom, w = A.readout(psi[:frame.k], frame)
            if not np.isfinite(C).all():
                return 1e3, np.zeros(P), {}
            r = C - C0
            f = 0.5 * float((r ** 2).sum())
            gC = r.ravel()
            d_dpsi = (frame.Wf @ gC - float(C.ravel() @ gC)) / denom
            J = 0.5 * circ.states_batch(th[None, :] + math.pi * np.eye(P))
            return f, J[:, :frame.k] @ d_dpsi, {}
        th, f, _, _ = A.adam(fg, th0, iters, lr)
        if best is None or f < best[0]:
            best = (f, th)
    th = best[1]
    C, denom, w = A.readout(circ.state(th)[:frame.k], frame)
    return th, C, rms(C - C0)


def grad_S_theta(circ, th, frame, sur):
    """d S~(C(theta)) / d theta, exact (state Jacobian by the pi shift, chain rule)."""
    P = th.size
    psi = circ.state(th)
    C, denom, w = A.readout(psi[:frame.k], frame)
    if not np.isfinite(C).all():
        return np.zeros(P)
    s, dS_dC = sur.value_grad(C)
    gC = dS_dC.ravel()
    d_dpsi = (frame.Wf @ gC - float(C.ravel() @ gC)) / denom
    J = 0.5 * circ.states_batch(th[None, :] + math.pi * np.eye(P))
    return J[:, :frame.k] @ d_dpsi


def circuit_step(circ, frame, sur, th, e, s_max=64.0):
    """One steepest-descent step on S~(C(theta)) in theta, scaled so the C-space RMS
    displacement from C(theta) is e (bisection on the scalar step).  NATIVE-FREE."""
    Cb, _, _ = A.readout(circ.state(th)[:frame.k], frame)
    g = grad_S_theta(circ, th, frame, sur)
    gn = np.linalg.norm(g)
    if gn < 1e-15:
        return Cb.copy(), 0.0, 0.0
    d = -g / gn

    def disp(s):
        C, _, _ = A.readout(circ.state(th + s * d)[:frame.k], frame)
        return rms(C - Cb) if np.isfinite(C).all() else float("inf")
    lo, hi = 0.0, 1e-3
    while disp(hi) < e and hi < s_max:
        lo, hi = hi, hi * 2.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if disp(mid) < e:
            lo = mid
        else:
            hi = mid
    s = 0.5 * (lo + hi)
    C, _, _ = A.readout(circ.state(th + s * d)[:frame.k], frame)
    return C, float(s), float(disp(s))


def run_ladder_target(pdb, seed=0):
    t0 = time.time()
    cand, frame, sur, C0, top = production_point(pdb)
    n = cand.n
    s0, g = sur.value_grad(C0)
    g = remove_rigid(g, C0)
    gr = rms(g)
    ghat = g / max(gr, 1e-15)
    clouds = {"prod": C0}
    row = dict(pdb=pdb, n=int(n), fold=int(cand.fold), fail18=bool(pdb in I.FAIL18), seed=int(seed),
               S_prod=float(s0), grad_rms=gr, arms={})

    def add(name, C, extra=None):
        r = dict(rmsd_cloud=oracle_rmsd(C, cand), S_smooth=float(sur.value_grad(C)[0]),
                 S_shipped=sur.shipped(C), disp_rms=rms(C - C0), **A.struct_diag(C))
        if extra:
            r.update(extra)
        row["arms"][name] = r
        clouds[name] = C
    add("prod", C0)
    for e in STEPS:
        add("step_e%g" % e, C0 - e * ghat)
        for d in range(N_RAND):
            v = SD.stable_rng(pdb, "s28A2_rand", int(d), int(seed), salt=SALT).normal(size=C0.shape)
            v = remove_rigid(v, C0); v = v / max(rms(v), 1e-15)
            add("rand%d_e%g" % (d, e), C0 + e * v)
    # the circuit family: nearest point to production, then one step
    circ = Q.StatevectorCircuit(A.N_QUBITS, A.LAYERS)
    thP, CP, resP = theta_nearest_production(circ, frame, C0)
    add("circP", CP, dict(residual_rms_to_prod=resP))
    for e in STEPS:
        C, s, achieved = circuit_step(circ, frame, sur, thP, e)
        add("circ_e%g" % e, C, dict(theta_step=s, disp_from_circP=achieved))
    row["secs"] = time.time() - t0
    A._save_structs(pdb, {k: v for k, v in clouds.items() if k != "prod"}, cand, suffix="_a2")
    return row


# ============================================================ driver
def phase(mode, pdbs, seed=0):
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"s28_A2_{mode}_rows.jsonl" if seed == 0 else f"s28_A2_{mode}_rows_seed{seed}.jsonl")
    done = A._jsonl_done(path)
    t0 = time.time()
    for k, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        t1 = time.time()
        if mode == "cosine":
            row = run_cosine_target(pdb)
            msg = " ".join("%s %+.3f" % (c, row["cos"][c]) for c in row["cos"]) + "  randref |cos| %.3f" % np.mean(np.abs(row["cos_random_ref"]))
        else:
            row = run_ladder_target(pdb, seed=seed)
            a = row["arms"]
            msg = "prod %.3f | step %.3f %.3f %.3f | rand0 %.3f %.3f %.3f | circP %.3f (res %.3f) circ %.3f %.3f %.3f" % (
                a["prod"]["rmsd_cloud"], a["step_e0.1"]["rmsd_cloud"], a["step_e0.3"]["rmsd_cloud"], a["step_e1"]["rmsd_cloud"],
                a["rand0_e0.1"]["rmsd_cloud"], a["rand0_e0.3"]["rmsd_cloud"], a["rand0_e1"]["rmsd_cloud"],
                a["circP"]["rmsd_cloud"], a["circP"]["residual_rms_to_prod"],
                a["circ_e0.1"]["rmsd_cloud"], a["circ_e0.3"]["rmsd_cloud"], a["circ_e1"]["rmsd_cloud"])
        A._append(path, row)
        print(f"  [{mode} {k+1}/{len(pdbs)}] {pdb} n={row['n']} {msg} {time.time()-t1:.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["cosine", "ladder"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    from s25 import phys_lib as P
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    phase(a.mode, pdbs, seed=a.seed)


if __name__ == "__main__":
    main()
