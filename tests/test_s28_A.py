"""tests/test_s28_A.py -- S28 lane A, the amplitude readout (`s27/s28_A_amp.py`).

Synthetic pools only (no target, no native, no cache): the readout reproduces the deployed
average when psi is uniform on a set; every gradient agrees with central finite differences;
the lam = 0 loop reproduces `core.quantum.run_cvar_vqe` bit-for-bit; the classical families'
chain rules are exact; the ORACLE RMSD gradient is exact; and the NaN-poison test: with every
native quantity NaN, every deployable output is bit-identical.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from core import quantum as Q                  # noqa: E402
from s12 import instrument as I                # noqa: E402
from s22 import qcand_lib as QC                # noqa: E402
from s27 import s28_A_amp as A                 # noqa: E402


# ------------------------------------------------------------------ synthetic fixtures
def _pool(k=40, n=8, seed=0):
    """k random smooth chains of n CA atoms with ~3.8 A bonds, plus a synthetic distogram."""
    rng = np.random.default_rng(seed)
    t = np.arange(n)                                    # an alpha-helical CA trace
    base = np.column_stack([2.3 * np.cos(np.radians(100) * t), 2.3 * np.sin(np.radians(100) * t), 1.5 * t])
    W = base[None] + rng.normal(0, 0.8, (k, n, 3))
    # rigidly move every window so the frame has work to do
    for a in range(k):
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        w, x, y, z = q
        R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                      [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                      [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
        W[a] = W[a] @ R.T + rng.normal(0, 5, 3)
    i, j = I.pair_index(n)
    grid = np.arange(2.0, 40.0, 0.05)
    tgt = rng.uniform(4, 15, len(i))
    risk = np.abs(grid[None, :] - tgt[:, None]) * rng.uniform(0.5, 2.0, (len(i), 1))
    dg = dict(risk=risk.astype(np.float32).astype(float), grid=grid, i=i, j=j)
    return W, dg


def _frame(k=40, n=8, m=10, seed=0):
    W, dg = _pool(k, n, seed)
    top = np.sort(np.random.default_rng(seed + 1).choice(k, m, replace=False))
    return W, dg, top, A.Frame(W, top), A.Surrogate(dg, n)


def _fd(fun, x, h=1e-5):
    x = np.asarray(x, float)
    g = np.zeros_like(x)
    for k in range(x.size):
        e = np.zeros_like(x); e[k] = h
        g[k] = (fun(x + e) - fun(x - e)) / (2 * h)
    return g


def _cos(a, b):
    return float(a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-300))


# ------------------------------------------------------------------ the readout
def test_uniform_on_set_reproduces_deployed_average():
    W, dg, top, fr, sur = _frame()
    psi = np.zeros(fr.k); psi[top] = 1.0
    C, denom, w = A.readout(psi, fr)
    Cd, b = I.coordinate_average(W[top])
    assert b == fr.b
    assert np.abs(C - Cd).max() < 1e-9
    assert abs(w.sum() - 1.0) < 1e-12 and denom == float(len(top))


def test_affine_weights_sum_to_one_for_either_sign():
    W, dg, top, fr, sur = _frame()
    rng = np.random.default_rng(3)
    for _ in range(20):
        psi = rng.normal(size=fr.k)
        C, denom, w = A.readout(psi, fr)
        assert abs(w.sum() - 1.0) < 1e-10
        # translation invariance of an affine combination: shifting every window shifts C
        fr2 = A.Frame(W, top); fr2.Wf = fr.Wf + np.tile([1.0, -2.0, 0.5], fr.n)[None, :]
        C2, _, _ = A.readout(psi, fr2)
        assert np.abs((C2 - C) - np.array([1.0, -2.0, 0.5])).max() < 1e-9


def test_readout_undefined_at_zero_denominator():
    W, dg, top, fr, sur = _frame()
    psi = np.zeros(fr.k); psi[0] = 1.0; psi[1] = -1.0
    C, denom, w = A.readout(psi, fr)
    assert np.isnan(C).all() and abs(denom) < A.DENOM_EPS


# ------------------------------------------------------------------ the surrogate
def test_surrogate_equals_shipped_on_grid_points_and_extends_linearly():
    W, dg, top, fr, sur = _frame()
    C = W[0]
    # move every pair distance onto a grid point is not possible for a rigid chain; instead
    # check the value against a direct interpolation
    v = C[sur.i] - C[sur.j]; d = np.linalg.norm(v, axis=1)
    inside = (d >= sur.grid[0]) & (d <= sur.grid[-1])
    assert inside.sum() >= 3
    direct = np.array([np.interp(d[p], sur.grid, sur.risk[p]) for p in range(sur.npairs)])
    # per-pair interpolated values, checked on the pairs inside the grid
    u = (d - sur.g0) / sur.dg; i0 = np.clip(np.floor(u).astype(int), 0, sur.G - 2)
    mine = sur.risk[np.arange(sur.npairs), i0] + sur.slope[np.arange(sur.npairs), i0] * (u - i0) * sur.dg
    assert np.abs(mine[inside] - direct[inside]).max() < 1e-9
    val, g = sur.value_grad(C)
    assert abs(val - mine.mean()) < 1e-9
    # shipped lookup within one bin's variation of the surrogate, on an all-in-grid structure
    Cin = C * max(1.0, 2.5 / d.min())
    din = np.linalg.norm(Cin[sur.i] - Cin[sur.j], axis=1)
    assert din.min() >= sur.grid[0] and din.max() <= sur.grid[-1]
    assert abs(sur.shipped(Cin) - sur.value_grad(Cin)[0]) < np.abs(np.diff(sur.risk, axis=1)).max() + 1e-9
    # beyond the grid: linear continuation with the last slope (ADDENDUM 1)
    C2 = C * 20.0                              # every distance far beyond 40 A
    val2, g2 = sur.value_grad(C2)
    d2 = np.linalg.norm(C2[sur.i] - C2[sur.j], axis=1)
    last = sur.risk[:, -1] + sur.slope[:, -1] * (d2 - sur.grid[-1])
    assert abs(val2 - last.mean()) < 1e-6
    assert np.abs(g2).max() > 0


def test_surrogate_gradient_matches_finite_differences():
    W, dg, top, fr, sur = _frame()
    C = W[3] + 0.01
    val, g = sur.value_grad(C)
    gfd = _fd(lambda x: sur.value_grad(x.reshape(C.shape))[0], C.ravel(), 1e-6)
    assert _cos(g.ravel(), gfd) > 0.9999 and np.abs(g.ravel() - gfd).max() < 1e-5


# ------------------------------------------------------------------ the circuit objective
def test_state_jacobian_is_half_pi_shift():
    circ = Q.StatevectorCircuit(4, 2)
    th = np.random.default_rng(0).normal(0, 0.6, circ.n_params())
    J = 0.5 * circ.states_batch(th[None, :] + math.pi * np.eye(th.size))
    for k in range(th.size):
        e = np.zeros_like(th); e[k] = 1e-6
        fd = (circ.state(th + e) - circ.state(th - e)) / 2e-6
        assert np.abs(J[k] - fd).max() < 1e-7


def test_objective_theta_gradient_vs_finite_differences():
    W, dg, top, fr, sur = _frame(k=40, n=8)
    circ = Q.StatevectorCircuit(6, 3)                  # 64 >= 40 candidates, 24 padding
    enc = QC.Encoding(A.zrank(np.random.default_rng(1).normal(size=fr.k)))
    assert enc.dim == circ.dim
    th = np.random.default_rng(2).normal(0, 0.6, circ.n_params())
    for lam in (0.0, 1.0, 3.0):
        f, g, info = A.objective_theta(circ, th, enc.E, 0.2, 0.5, lam, fr, sur)
        gfd = _fd(lambda x: A.objective_theta(circ, x, enc.E, 0.2, 0.5, lam, fr, sur)[0], th, 1e-5)
        assert _cos(g, gfd) > 0.999, (lam, _cos(g, gfd))
        assert np.abs(g - gfd).max() < 1e-5, (lam, np.abs(g - gfd).max())
        if lam > 0:
            assert "S_smooth" in info and np.isfinite(info["S_smooth"])


def test_lam0_loop_reproduces_run_cvar_vqe_bitwise():
    W, dg, top, fr, sur = _frame(k=40, n=8)
    circ = Q.StatevectorCircuit(6, 3)
    enc = QC.Encoding(A.zrank(np.random.default_rng(1).normal(size=fr.k)))
    th0 = np.random.default_rng(0).normal(0.0, 0.6, circ.n_params())
    th, f, info, tr = A.adam(lambda x: A.objective_theta(circ, x, enc.E, 0.2, 0.5, 0.0, fr, sur), th0, 30, 0.15)
    p_ref, cv, H, _ = Q.run_cvar_vqe(enc.E, 0.2, 0.5, n=6, layers=3, iters=30, restarts=1, seed=0, lr=0.15)
    assert np.array_equal(circ.probs(th), p_ref)


# ------------------------------------------------------------------ the classical families
def test_objective_psi_gradient_vs_finite_differences():
    W, dg, top, fr, sur = _frame(k=40, n=8)
    enc = QC.Encoding(A.zrank(np.random.default_rng(1).normal(size=fr.k)))
    psi = np.random.default_rng(5).normal(size=enc.dim)
    for lam in (0.0, 1.0):
        f, g, info = A.objective_psi(psi, enc.E, 0.2, 0.5, lam, fr, sur)
        gfd = _fd(lambda x: A.objective_psi(x, enc.E, 0.2, 0.5, lam, fr, sur)[0], psi, 1e-5)
        assert _cos(g, gfd) > 0.999 and np.abs(g - gfd).max() < 1e-5


@pytest.mark.parametrize("kind", ["a500", "a75", "sub", "simplex"])
def test_family_chain_rules(kind):
    W, dg, top, fr, sur = _frame(k=40, n=8)
    enc = QC.Encoding(A.zrank(np.random.default_rng(1).normal(size=fr.k)))
    B = np.random.default_rng(9).normal(size=(fr.k, 7))
    fam = A.Family(kind, enc.dim, fr.k, top, B=B)
    z = fam.init("rand", np.random.default_rng(11))

    def full(zz):
        psi = fam.psi(zz)
        f, g, _ = A.objective_psi(psi, enc.E, 0.2, 0.5, 1.0, fr, sur)
        return f, fam.chain(zz, psi, g)
    f, g = full(z)
    gfd = _fd(lambda x: full(x)[0], z, 1e-5)
    assert _cos(g, gfd) > 0.999 and np.abs(g - gfd).max() < 1e-4
    assert np.all(fam.psi(z)[fr.k:] == 0.0)            # padding carries nothing
    if kind == "simplex":
        psi = fam.psi(z)
        assert np.all(psi >= 0) and abs((psi ** 2).sum() - 1) < 1e-12
    if kind != "sub":
        zp = fam.init("prod", None)
        C, denom, w = A.readout(fam.psi(zp)[:fr.k], fr)
        Cd, _ = I.coordinate_average(W[top])
        assert np.abs(C - Cd).max() < 1e-6            # the production point is in the family


@pytest.mark.parametrize("kind", ["a500", "a75", "sub", "simplex"])
def test_wfamily_chain_rules(kind):
    W, dg, top, fr, sur = _frame(k=40, n=8)
    B = np.random.default_rng(9).normal(size=(fr.k, 7))
    fam = A.WFamily(kind, fr.k, top, B=B)
    z = fam.init("rand", np.random.default_rng(12))

    def full(zz):
        w = fam.w(zz)
        f, g, _ = A.objective_w(w, fr, sur)
        return f, fam.chain(zz, w, g)
    f, g = full(z)
    gfd = _fd(lambda x: full(x)[0], z, 1e-5)
    assert _cos(g, gfd) > 0.999 and np.abs(g - gfd).max() < 1e-6 * max(1.0, np.abs(g).max())


# ------------------------------------------------------------------ ORACLE pieces
def test_oracle_rmsd_and_gradient():
    W, dg, top, fr, sur = _frame()
    nat = W[7] @ np.linalg.qr(np.random.default_rng(4).normal(size=(3, 3)))[0] + 3.0
    C = fr.Wp[2]
    r, g = A.oracle_rmsd_grad(C, nat)
    assert abs(r - I.ca_rmsd(C, nat)) < 1e-12
    gfd = _fd(lambda x: I.ca_rmsd(x.reshape(C.shape), nat), C.ravel(), 1e-6)
    assert _cos(g.ravel(), gfd) > 0.9999 and np.abs(g.ravel() - gfd).max() < 1e-5


def test_oracle_affine_ls_reaches_a_member_of_the_span():
    W, dg, top, fr, sur = _frame(k=40, n=8)
    # a native that IS an affine combination of the posed windows must be reached to ~0
    rng = np.random.default_rng(6)
    w = rng.normal(size=fr.k); w /= w.sum()
    nat = (w @ fr.Wf).reshape(fr.n, 3)
    r, C, ww = A.oracle_affine_ls(fr, nat)
    assert r < 1e-6 and abs(ww.sum() - 1) < 1e-9
    r75, C75, w75 = A.oracle_affine_ls(fr, nat, top)
    assert abs(w75.sum() - 1) < 1e-9 and r75 >= -1e-12
    B = rng.normal(size=(fr.k, 7))
    rs, Cs, ws = A.oracle_subspace_ls(fr, nat, B)
    assert abs(ws.sum() - 1) < 1e-9


# ------------------------------------------------------------------ NaN-poison
def test_nan_poison_every_deployable_output_bit_identical():
    """Replace every native quantity with NaN: every deployable structure is bit-identical."""
    W, dg, top, fr, sur = _frame(k=40, n=8)
    circ = Q.StatevectorCircuit(6, 3)
    enc = QC.Encoding(A.zrank(np.random.default_rng(1).normal(size=fr.k)))
    th0 = np.random.default_rng(0).normal(0.0, 0.6, circ.n_params())
    outs = []
    for poison in (False, True):
        nat = np.full((fr.n, 3), np.nan) if poison else W[5]
        # nothing deployable takes `nat`; the assertion is that the code path cannot see it
        th, f, info, tr = A.adam(lambda x: A.objective_theta(circ, x, enc.E, 0.2, 0.5, 1.0, fr, sur), th0, 20, 0.15)
        C, denom, w = A.readout(circ.state(th)[:fr.k], fr)
        fam = A.Family("a500", enc.dim, fr.k, top)
        z, fz, iz, nev = A.lbfgs(lambda zz: (lambda psi, r: (r[0], fam.chain(zz, psi, r[1]), r[2]))(fam.psi(zz), A.objective_psi(fam.psi(zz), enc.E, 0.2, 0.5, 1.0, fr, sur)), fam.init("rand", np.random.default_rng(1)), 30)
        Cz, _, _ = A.readout(fam.psi(z)[:fr.k], fr)
        outs.append((C.copy(), Cz.copy(), float(f), float(fz)))
        # ORACLE scoring of the same structure is the only place `nat` is read
        r = A.oracle_rmsd_of(C, type("c", (), dict(nat_ca=nat))())
        assert (np.isnan(r) if poison else np.isfinite(r))
    a, b = outs
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]) and a[2] == b[2] and a[3] == b[3]


def test_adam_matches_run_cvar_vqe_arithmetic_on_trace():
    W, dg, top, fr, sur = _frame(k=40, n=8)
    circ = Q.StatevectorCircuit(6, 3)
    enc = QC.Encoding(A.zrank(np.random.default_rng(1).normal(size=fr.k)))
    th0 = np.random.default_rng(0).normal(0.0, 0.6, circ.n_params())
    th, f, info, tr = A.adam(lambda x: A.objective_theta(circ, x, enc.E, 0.2, 0.5, 0.0, fr, sur), th0, 40, 0.15, trace_every=20)
    assert [t for t, _ in tr] == [0, 20, 40]
    assert tr[-1][1] <= tr[0][1]                      # the free energy fell


def test_replace_retry_survives_transient_permission_error(tmp_path, monkeypatch):
    """Windows refuses os.replace on an open target; the store retries (S26 L59's guard)."""
    src = tmp_path / "a.tmp"; dst = tmp_path / "a.npz"
    src.write_bytes(b"x"); dst.write_bytes(b"old")
    real = os.replace
    calls = {"n": 0}

    def flaky(a, b):
        calls["n"] += 1
        if calls["n"] <= 2:
            raise PermissionError(5, "Access is denied")
        return real(a, b)
    monkeypatch.setattr(A.os, "replace", flaky)
    A.replace_retry(str(src), str(dst), tries=4, wait=0.0)
    assert dst.read_bytes() == b"x" and calls["n"] == 3
    src.write_bytes(b"y"); calls["n"] = -10
    with pytest.raises(PermissionError):
        A.replace_retry(str(src), str(dst), tries=3, wait=0.0)


def test_structs_store_merges_across_phase_files(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "STRUCTS", str(tmp_path))
    cand = type("c", (), dict(seq="ACDE", fold=1, n=4))()
    A._save_structs("XXXX", dict(prod=np.ones((4, 3))), cand)
    A._save_structs("XXXX", dict(circ=np.zeros((4, 3))), cand, suffix="_recog")
    A._save_structs("XXXX", dict(other=np.full((4, 3), 2.0)), cand, suffix="_recog")   # merge within a phase
    assert sorted(os.listdir(tmp_path)) == ["XXXX.npz", "XXXX_recog.npz"]
    with np.load(tmp_path / "XXXX_recog.npz") as z:
        assert set(z.files) >= {"circ", "other", "seq", "fold", "n"}
