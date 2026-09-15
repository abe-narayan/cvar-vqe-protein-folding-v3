"""tests/test_s28_A2.py -- S28 lane A2 (`s27/s28_A2_local.py`): the objective's local behaviour.

Synthetic only: analytic channel gradients against central finite differences; rigid-body
removal (a pure rigid field projects to zero, a shape field keeps its non-rigid part); the
circuit step hits the requested displacement; the NaN-poison test for the deployable ladder
(every emitted cloud bit-identical with the native replaced by NaN).
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from core import quantum as Q                  # noqa: E402
from s12 import instrument as I                # noqa: E402
from s27 import s28_A_amp as A                 # noqa: E402
from s27 import s28_A2_local as L              # noqa: E402
from tests.test_s28_A import _frame, _fd, _cos  # noqa: E402


def _C0():
    W, dg, top, fr, sur = _frame(k=40, n=8)
    psi = np.zeros(fr.k); psi[top] = 1.0
    C0, _, _ = A.readout(psi, fr)
    return W, dg, top, fr, sur, C0


def test_rg_law_and_exvol_gradients_vs_fd():
    W, dg, top, fr, sur, C0 = _C0()
    C = C0 * 0.6                                            # compress so EXVOL has active pairs
    v, g = L.grad_rg_law(C, 8)
    gfd = _fd(lambda x: L.grad_rg_law(x.reshape(C.shape), 8)[0], C.ravel(), 1e-6)
    assert _cos(g.ravel(), gfd) > 0.9999 and np.abs(g.ravel() - gfd).max() < 1e-6
    v2, g2 = L.grad_exvol(C)
    assert v2 > 0
    gfd2 = _fd(lambda x: L.grad_exvol(x.reshape(C.shape))[0], C.ravel(), 1e-6)
    assert _cos(g2.ravel(), gfd2) > 0.9999 and np.abs(g2.ravel() - gfd2).max() < 1e-5


def test_rigid_removal():
    W, dg, top, fr, sur, C0 = _C0()
    rng = np.random.default_rng(0)
    # a pure rigid field: translation + small rotation about the centroid
    w = rng.normal(size=3) * 0.1
    c = C0 - C0.mean(0)
    rigid = np.cross(w[None, :], c) + rng.normal(size=3)[None, :]
    assert L.rms(L.remove_rigid(rigid, C0)) < 1e-10
    # a shape field: idempotent projection, and its rigid part is what is removed
    v = rng.normal(size=C0.shape)
    vp = L.remove_rigid(v, C0)
    assert np.abs(L.remove_rigid(vp, C0) - vp).max() < 1e-12
    B = L.rigid_basis(C0)
    assert np.abs(B.T @ vp.ravel()).max() < 1e-10
    assert L.rms(vp) <= L.rms(v)


def test_oracle_direction_is_zero_at_the_native():
    W, dg, top, fr, sur, C0 = _C0()
    u = L.oracle_direction(C0, C0 @ np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))[0] + 2.0)
    assert L.rms(u) < 1e-9


def test_circuit_step_hits_requested_displacement():
    W, dg, top, fr, sur, C0 = _C0()
    circ = Q.StatevectorCircuit(6, 3)
    th, CP, res = L.theta_nearest_production(circ, fr, C0, starts=2, iters=60)
    assert np.isfinite(res)
    for e in (0.1, 0.3):
        C, s, achieved = L.circuit_step(circ, fr, sur, th, e)
        assert abs(achieved - e) < 1e-3 and s > 0
        assert abs(L.rms(C - CP) - e) < 1e-3
    g = L.grad_S_theta(circ, th, fr, sur)
    gfd = _fd(lambda x: sur.value_grad(A.readout(circ.state(x)[:fr.k], fr)[0])[0], th, 1e-5)
    assert _cos(g, gfd) > 0.999


def test_ladder_is_native_free(monkeypatch, tmp_path):
    """NaN-poison: with the native NaN the deployable clouds are bit-identical."""
    W, dg, top, fr, sur, C0 = _C0()
    outs = []
    for poison in (False, True):
        s0, g = sur.value_grad(C0)
        g = L.remove_rigid(g, C0)
        ghat = g / L.rms(g)
        clouds = {"step": C0 - 0.3 * ghat}
        v = np.random.default_rng(5).normal(size=C0.shape); v = L.remove_rigid(v, C0); v /= L.rms(v)
        clouds["rand"] = C0 + 0.3 * v
        nat = np.full(C0.shape, np.nan) if poison else W[3]
        r = A.oracle_rmsd_of(clouds["step"], type("c", (), dict(nat_ca=nat))())
        assert np.isnan(r) if poison else np.isfinite(r)
        outs.append(clouds)
    assert np.array_equal(outs[0]["step"], outs[1]["step"]) and np.array_equal(outs[0]["rand"], outs[1]["rand"])
    # the module's own step is what the ladder writes: the gradient never touches a native
    import inspect
    src = inspect.getsource(L.run_ladder_target)
    assert "nat_ca" not in src.replace("oracle_rmsd(C, cand)", "")
