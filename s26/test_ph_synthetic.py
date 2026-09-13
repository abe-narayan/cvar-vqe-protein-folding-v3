"""s26/test_ph_synthetic.py -- unit tests for the PH lane on SYNTHETIC data only.

No native, no RMSD to a native, no artefact read.  Run with
    python s26/test_ph_synthetic.py
"""
from __future__ import annotations

import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s26 import ph_lib as L                                   # noqa: E402
from s26 import ph_reject as R                                # noqa: E402
from s26 import ph_c3 as C3                                   # noqa: E402


def test_gate():
    with tempfile.TemporaryDirectory() as d:
        closed = os.path.join(d, "closed.md"); open_ = os.path.join(d, "open.md")
        with open(closed, "w") as fh:
            # the decoy: ledger L5 quotes the phrase inside a sentence and must NOT open the gate
            fh.write("## L1 -- nothing\nThe gate forbids any endpoint experiment before\n"
                     "\"PHASE 0 SIGNED OFF\". The census is allowed.\n")
        with open(open_, "w") as fh:
            fh.write("## L9 -- PHASE 0 SIGNED OFF (date, coordinator)\n")
        assert not L.gate_open(closed) and L.gate_open(open_)
        with open(open_, "w") as fh:
            fh.write("## L9 -- gate\n\nPHASE 0 SIGNED OFF\n")
        assert L.gate_open(open_)
        try:
            L.require_gate("x", ledger=closed)
            raise AssertionError("gate did not raise")
        except RuntimeError:
            pass
        L.require_gate("x", ledger=open_)
    # oracle keys are refused by the native-free universe loader
    try:
        L.univ_nativefree("XXXX", keys=("W", "rr"))
        raise AssertionError("oracle key accepted")
    except ValueError:
        pass


def test_omega():
    from core import geometry as geo
    n = 12
    phi = np.full(n, math.radians(-63.0)); psi = np.full(n, math.radians(-42.0))
    bb = geo.build_backbone(phi, psi)
    om = L.omega_deg(bb["CA"], bb["C"], bb["N"])
    assert om.shape == (n - 1,) and np.allclose(np.abs(om), 180.0, atol=1e-6), om
    st = L.consecutive_ca(bb["CA"])
    assert np.allclose(st, 3.804, atol=2e-3), st
    cis = geo.build_backbone(phi, psi, omega=0.0)
    omc = L.omega_deg(cis["CA"], cis["C"], cis["N"])
    stc = L.consecutive_ca(cis["CA"])
    assert np.allclose(np.abs(omc), 0.0, atol=1e-6)
    assert (stc < L.CIS_CA_CA).all() and (stc > 2.7).all(), stc
    # the two census criteria agree on ideal cis and ideal trans
    assert ((np.abs(omc) < L.CIS_OMEGA_DEG) == (stc < L.CIS_CA_CA)).all()
    assert ((np.abs(om) < L.CIS_OMEGA_DEG) == (st < L.CIS_CA_CA)).all()


def test_reject_operators():
    rng = np.random.default_rng(0)
    K, m = 500, 75
    e = 10 ** rng.uniform(-3, 9, K)
    score = rng.normal(size=K)
    order = L.argsort_stable(score)
    sub = order[:m]
    T = 1e4
    S = R.reject_shrink(sub, e, T)
    assert (e[S] <= T).all() and set(S) <= set(sub)
    Rt = R.reject_refill(order, e, T)
    # R equals the top-m of the survivors in score order (pre-filter form)
    surv = order[e[order] <= T][:m]
    assert np.array_equal(Rt, surv) and len(Rt) == m
    # every retained-but-not-in-sub member ranks below every sub member that survived
    rank = {int(i): k for k, i in enumerate(order)}
    assert max(rank[int(i)] for i in Rt) >= m - 1
    r = m - len(S)
    for _ in range(5):
        rs = R.random_shrink(sub, r, rng)
        assert len(rs) == len(S) and set(rs) <= set(sub)
        rr = R.random_refill(order, sub, r, rng)
        assert len(rr) == m and len(set(rr.tolist())) == m
        new = set(rr.tolist()) - set(sub.tolist())
        assert all(rank[i] >= m for i in new) and len(new) == r
    ep = R.permuted_energy(e, rng)
    assert np.array_equal(np.sort(ep), np.sort(e)) and not np.array_equal(ep, e)
    # a threshold above every energy rejects nothing and R == anchor as a set
    R0 = R.reject_refill(order, e, 1e12)
    assert set(R0.tolist()) == set(sub.tolist())
    # a threshold below every energy: S empty, R empty
    assert len(R.reject_shrink(sub, e, 1e-9)) == 0 and len(R.reject_refill(order, e, 1e-9)) == 0
    pool = {"pdb": "TEST", "e": e, "order": order, "sub": sub}
    rs = R.retained_sets(pool, T, n_draws=4)
    assert rs["n_reject"] == r and len(rs["RANDS"]) == 4 and len(rs["PERMR"]) == 4
    assert all(len(x) == len(S) for x in rs["RANDS"]) and all(len(x) == m for x in rs["RANDR"])


def test_random_displacement():
    rng = np.random.default_rng(1)
    n = 14
    ca = np.cumsum(rng.normal(size=(n, 3)), 0) * 2.0
    mag = 0.37
    g = L.random_displacement(ca, mag, rng)
    assert abs(math.sqrt((g ** 2).sum(1).mean()) - mag) < 1e-9
    Q = L.rigid_basis(ca)
    assert np.abs(Q.T @ g.ravel()).max() < 1e-9, "rigid components not removed"
    # the identity of S16 section 3.4 is an upper bound on the superposed RMSD
    nat = ca + rng.normal(size=(n, 3)) * 1.5
    out = ca + g
    v, m2 = C3.displacement(ca, out)
    assert abs(m2 - mag) < 1e-6
    r = L.superpose_onto(nat, ca) - ca
    pred = math.sqrt(max(((r ** 2).sum() - 2 * (v * r).sum() + (v ** 2).sum()) / n, 0.0))
    real = L.ca_rmsd(out, nat)
    assert real <= pred + 1e-9, (real, pred)


def test_ca_helpers():
    ca = np.zeros((10, 3)); ca[:, 0] = np.arange(10) * 3.8
    assert L.ca_contacts(ca, 4.0) == 0 and L.min_sep3(ca) > 11.0
    assert abs(L.rg_of(ca) - math.sqrt(((np.arange(10) - 4.5) ** 2).mean()) * 3.8) < 1e-9


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("ALL PH SYNTHETIC TESTS PASS")
