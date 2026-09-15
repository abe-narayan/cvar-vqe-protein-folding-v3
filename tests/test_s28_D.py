"""tests/test_s28_D.py -- lane D (the Adversary): the checks the lanes' own test files do not make.

1. Lane A's `test_nan_poison_every_deployable_output_bit_identical` never passes a poisoned pool
   through the deployable path (the `nat` it builds is read only by `oracle_rmsd_of`); it is a
   determinism test.  Here the WHOLE recognition phase (`run_recog_target`) runs twice on a
   synthetic pool, once with the natives intact and once with `nat_ca` / `oracle_rr` NaN, through
   a monkeypatched `load_pool`, with the run constants shrunk so it takes seconds; every emitted
   structure, every objective value and every weight vector must be bit-identical, and every
   ORACLE `rmsd_cloud` must be NaN under poison.
2. The same end to end through the built-chain phase (`run_chain_target`) with a stubbed
   projection: the chain coordinates are identical and the chain RMSD is NaN under poison.
3. Lane B's sign coherence and Perron reading: the exact ground state of diag(E) - J A with
   A >= 0 entrywise and connected has one sign (so `sign_coherence` is 1.0), and `hop_value`
   at the Perron vector of a unit-spectral-norm A is 1.0 (the bound the check in S28-L2 uses).
4. Lane C's readouts: the matched control changes the answer for k < 74 / q > 0 / beta > 0 on
   a pool where the ranker is informative, and the identity cells do not (so a "control equals
   arm" reading cannot arise from a no-op control).

Synthetic data only; nothing on disk is read except the frozen modules.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from core import quantum as Q              # noqa: E402
from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s27 import s28_A_amp as A             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402
from s27 import s28_C_readout as R         # noqa: E402


# ------------------------------------------------------------------ synthetic pool + distogram
def _synthetic_target(k=40, n=8, seed=0):
    rng = np.random.default_rng(seed)
    base = np.cumsum(rng.normal(0.0, 1.0, (n, 3)), axis=0) * 2.0
    W = base[None] + rng.normal(0.0, 1.0, (k, n, 3))
    nat = base + rng.normal(0.0, 0.5, (n, 3))
    i, j = I.pair_index(n)
    grid = np.arange(2.0, 22.0, 0.05)
    G = len(grid)
    risk = np.abs(grid[None, :] - rng.uniform(4, 12, (len(i), 1))) * rng.uniform(0.5, 1.5, (len(i), 1))
    dg = {"risk": risk, "grid": grid, "i": i, "j": j}
    dis = rng.normal(size=k)
    return W, nat, dg, dis


class _Cand:
    def __init__(self, W, nat, rr, pdb="SYND", seq="ACDEFGHI", fold=0):
        self.W = np.asarray(W, float); self.nat_ca = nat; self.oracle_rr = rr
        self.pdb, self.seq, self.fold = pdb, seq, int(fold)
        self.k, self.n = int(self.W.shape[0]), int(self.W.shape[1])


def _shrink(monkeypatch):
    monkeypatch.setattr(A, "N_QUBITS", 6)
    monkeypatch.setattr(A, "ITERS", 6)
    monkeypatch.setattr(A, "ITERS_LONG", 6)
    monkeypatch.setattr(A, "LAMS", (1.0,))
    monkeypatch.setattr(A, "N_SUB", 1)
    monkeypatch.setattr(A, "N_UNTRAINED", 2)
    monkeypatch.setattr(A, "BUDGET_MATCHED", 6)
    monkeypatch.setattr(A, "BUDGET_CONV", 8)
    monkeypatch.setattr(A, "_save_structs", lambda *a, **k: None)


def _patched_pool(monkeypatch, poison):
    W, nat, dg, dis = _synthetic_target()
    k = len(W)
    rr = np.linalg.norm(W - nat[None], axis=(1, 2))
    if poison:
        nat = np.full_like(nat, np.nan)
        rr = np.full(k, np.nan)
    cand = _Cand(W, nat, rr)
    top = np.sort(np.argsort(dis, kind="stable")[:15])
    monkeypatch.setattr(A, "load_pool", lambda pdb: (cand, dis.copy(), top, dg))
    monkeypatch.setattr(A, "M", 15)
    return cand


def _flatten(row):
    out = {}
    for name, r in row["arms"].items():
        for key, v in r.items():
            if key in ("secs",):
                continue
            out[(name, key)] = v
    return out


def test_lane_A_recognition_phase_is_bit_identical_under_nan_poison(monkeypatch):
    """The real poison test for lane A: the deployable path never sees the native."""
    _shrink(monkeypatch)
    rows = []
    for poison in (False, True):
        _patched_pool(monkeypatch, poison)
        rows.append(A.run_recog_target("SYND", seed=0, long_diag=True))
    a, b = _flatten(rows[0]), _flatten(rows[1])
    assert set(a) == set(b)
    n_checked = 0
    for key in a:
        name, field = key
        va, vb = a[key], b[key]
        if field == "rmsd_cloud":
            assert np.isfinite(va) or rows[0]["arms"][name]["undefined"], key
            assert np.isnan(vb), ("poison did not reach the ORACLE score", key)
            continue
        if isinstance(va, (list, tuple, np.ndarray)):
            assert np.array_equal(np.asarray(va, float), np.asarray(vb, float), equal_nan=True), key
        elif isinstance(va, float):
            assert (va == vb) or (np.isnan(va) and np.isnan(vb)), key
        else:
            assert va == vb, key
        n_checked += 1
    assert n_checked > 50
    assert "circ_l1_i6" in rows[0]["arms"] and "untr_0" in rows[0]["arms"]


def test_lane_A_chain_phase_is_bit_identical_under_nan_poison(monkeypatch, tmp_path):
    """The built-chain phase: identical projected coordinates under poison.

    `run_chain_target` scores the chain with `I.ca_rmsd(ca, nat)` in the same dict it builds
    the chain in, and `kabsch_rmsd_batch` RAISES (SVD) on a NaN native rather than returning
    NaN, so the phase cannot be run to completion under poison as written.  The deployable
    half is captured through a shim on the ORACLE scorer only: it records every `ca` it is
    handed and returns NaN when the native is non-finite.  Nothing deployable is patched.
    """
    W, nat, dg, dis = _synthetic_target()
    n = W.shape[1]
    f = tmp_path / "SYND.npz"
    np.savez_compressed(f, prod=W[:10].mean(0), circ=W[3], seq=np.array("ACDEFGHI"), fold=np.array(0), n=np.array(n))
    monkeypatch.setattr(A, "STRUCTS", str(tmp_path))
    monkeypatch.setattr(I, "project", lambda C, seq, fold: {"ca": np.asarray(C, float) * 1.01 + 0.1,
                                                            "phi": np.zeros(n), "psi": np.zeros(n), "fit_ca": np.asarray(C, float)})
    real_rmsd = I.ca_rmsd
    seen = []

    def shim(a, b):
        seen.append(np.asarray(a, float).copy())
        if not np.isfinite(np.asarray(b, float)).all():
            return float("nan")
        return real_rmsd(a, b)

    monkeypatch.setattr(I, "ca_rmsd", shim)
    outs = []
    for poison in (False, True):
        u = {"nat_ca": np.full((n, 3), np.nan) if poison else nat}
        monkeypatch.setattr(I, "load_univ", lambda pdb, u=u: u)
        outs.append(A.run_chain_target("SYND", ["prod", "circ"]))
    a, b = outs[0]["arms"], outs[1]["arms"]
    for arm in ("prod", "circ"):
        assert np.isfinite(a[arm]["rmsd_chain"]) and np.isnan(b[arm]["rmsd_chain"])
        assert a[arm]["rg"] == b[arm]["rg"] and a[arm]["bond"] == b[arm]["bond"]
    # the scorer saw the same chain coordinates in both runs (2 arms x 2 calls each)
    assert len(seen) == 8
    for i in range(4):
        assert np.array_equal(seen[i], seen[4 + i])


def test_lane_B_perron_state_is_sign_coherent_and_attains_the_hop_bound():
    rng = np.random.default_rng(5)
    k = 32
    X = rng.normal(size=(k, 3))
    D = np.linalg.norm(X[:, None] - X[None], axis=-1)
    g = B.kernel_graph(D)
    Aa = g["A"]
    assert np.allclose(Aa, Aa.T) and np.all(np.diag(Aa) == 0) and np.all(Aa >= 0)
    w, v = np.linalg.eigh(Aa)
    perron = v[:, -1] * np.sign(v[:, -1].sum())
    assert abs(w[-1] - 1.0) < 1e-9                       # unit spectral norm
    assert abs(B.hop_value(perron, Aa) - 1.0) < 1e-9     # the bound S28-L2(b) quotes
    assert abs(B.sign_coherence(perron) - 1.0) < 1e-9    # one sign
    E = B.deployed_E(k)[rng.permutation(k)]
    gs = B.ground_state(E, Aa, 3.0)
    assert B.sign_coherence(gs["psi"]) > 0.999           # Perron-like at large J
    assert np.all(gs["psi"] > -1e-12)
    # a random RY state does not have one sign: the circuit has to find it
    circ = Q.StatevectorCircuit(5, 3)
    th = np.random.default_rng(1).normal(0.0, 0.6, circ.n_params())
    assert B.sign_coherence(circ.state(th)) < 0.9


def test_lane_C_matched_controls_are_not_no_ops():
    rng = np.random.default_rng(9)
    k, n = 80, 10
    base = np.cumsum(rng.normal(size=(n, 3)), 0) * 2
    W = base[None] + rng.normal(0.0, 1.5, (k, n, 3))
    W[:20] += rng.normal(0.0, 2.0, (20, 1, 3))          # a distinct cluster
    dis = rng.normal(size=k)
    cons = np.linalg.norm(W - W.mean(0)[None], axis=(1, 2)) + 0.01 * rng.normal(size=k)
    ch = {"DIS": dis, "CONS": cons, "DISTPOT": rng.normal(size=k)}
    cand = _Cand(W, None, None, pdb="SYNC")
    import s27.run_pool as RP
    top = RP.topm(dis, R.M, RP.rng_for("SYNC", "tiekey").random(k))
    prod, _ = R.emit(cand, ch, None, top=top)
    for real, perm, identity in ((("mednb", "CONS", False, 20), ("mednb", "CONS", True, 20), ("mednb", "CONS", False, 74)),
                                 (("trim", "CONS", False, 0.1), ("trim", "CONS", True, 0.1), ("trim", "CONS", False, 0.0)),
                                 (("divw", "CONS", False, (1.0, 1.0)), ("divw", "CONS", True, (1.0, 1.0)), ("divw", "CONS", False, (0.0, 0.0)))):
        Cr, _ = R.emit(cand, ch, real, top=top)
        Cp, _ = R.emit(cand, ch, perm, top=top)
        Ci, _ = R.emit(cand, ch, identity, top=top)
        assert not np.allclose(Cr, Cp), real            # the control moves the answer
        assert np.abs(Ci - prod).max() < 1e-9, identity  # the identity cell is production


def test_lane_C_trim_and_mednb_never_read_array_order_for_ties():
    """A ranker with EXACT ties: the retained set depends on the stable key, not on the order
    of the ties in the array (permuting tied entries with the key fixed changes nothing)."""
    rng = np.random.default_rng(2)
    k, n = 75, 8
    W = rng.normal(size=(k, n, 3))
    P = I.pairwise_rmsd(W)
    r_tied = np.repeat(np.arange(15.0), 5)              # 15 groups of 5 exact ties
    import s27.run_pool as RP
    key = RP.rng_for("SYNT", "s28C-tie-CONS").random(k)
    order = np.lexsort((key, r_tied))
    rr = np.empty(k); rr[order] = np.arange(k, dtype=float)
    _, keep1 = R.readout_trim(W, P, rr, 0.1)
    # the same ties presented in reversed array order, the same key per member
    perm = np.arange(k)[::-1]
    order2 = np.lexsort((key[perm], r_tied[perm]))
    rr2 = np.empty(k); rr2[order2] = np.arange(k, dtype=float)
    _, keep2 = R.readout_trim(W[perm], P[np.ix_(perm, perm)], rr2, 0.1)
    assert set(perm[keep2].tolist()) == set(keep1.tolist())
