"""Tests for `s27/s28_B_hop.py` (S28 lane B: H = diag(E) - J A inside the CVaR-VQE).

What each test pins:
  * the hopping gradient by parameter shift equals central finite differences (independent
    machinery), on a random symmetric A at n = 5 and on a padded graph at n = 6;
  * the graph is symmetric with zero diagonal, unit spectral norm, zero padding rows/cols;
  * the PERM control keeps the spectrum and the degree multiset; the RAND control keeps the
    degree SEQUENCE (row sums to 1e-8) and is symmetric with zero diagonal;
  * J = 0 reproduces `core.quantum.run_cvar_vqe` bit for bit (`==` on p), hence
    `s24.d_harness.arm_vqe`; on one real target when the S27 cache is present;
  * the exact ground state is the argmin at J = 0 and satisfies the eigen-equation at J > 0;
  * for J > 0 the objective is not a function of p: two states with equal p and different
    signs give different F;
  * NaN-poison: the native-free half's output is bit-identical when nat_ca / oracle_rr are NaN;
  * the R3 readout breaks ties by the key, not by array order.

Run:  python -m pytest tests/test_s28_B.py -q -p no:cacheprovider
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from core import quantum as Q                      # noqa: E402
from s27 import s28_B_hop as B                     # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HAVE_CACHE = os.path.exists(os.path.join(ROOT, "s27", "cache", "1A13.npz")) and \
    os.path.exists(os.path.join(ROOT, "s24", "cache_amber", "1A13.npz"))


def _rand_sym(dim, rng):
    A = rng.random((dim, dim))
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def _cos(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


# ------------------------------------------------------------------ the hopping gradient
def test_hop_gradient_paramshift_matches_finite_differences_n5():
    rng = np.random.default_rng(3)
    circ = Q.StatevectorCircuit(5, 3)
    A = _rand_sym(circ.dim, rng)
    for _ in range(3):
        th = rng.normal(0.0, 0.6, circ.n_params())
        g_ps = B.grad_hop_paramshift(circ, th, A)
        g_fd = B.grad_hop_fd(circ, th, A)
        assert g_ps.shape == (circ.n_params(),)
        assert np.max(np.abs(g_ps - g_fd)) < 1e-6
        assert _cos(g_ps, g_fd) > 1 - 1e-9


def test_hop_gradient_on_padded_graph_n6():
    rng = np.random.default_rng(7)
    circ = Q.StatevectorCircuit(6, 3)
    W = rng.normal(0.0, 3.0, (50, 8, 3))
    D = B.pairwise_rmsd_matrix(W)
    g = B.kernel_graph(D)
    A = B.pad_graph(g["A"], circ.dim)
    th = rng.normal(0.0, 0.6, circ.n_params())
    g_ps = B.grad_hop_paramshift(circ, th, A)
    g_fd = B.grad_hop_fd(circ, th, A)
    assert np.max(np.abs(g_ps - g_fd)) < 1e-6


def test_full_objective_gradient_matches_finite_differences():
    """The whole F = CVaR - T H - J <A>, differentiated two ways."""
    rng = np.random.default_rng(11)
    circ = Q.StatevectorCircuit(5, 3)
    E = B.deployed_E(circ.dim)
    A = _rand_sym(circ.dim, rng)
    A = A / np.linalg.eigvalsh(A)[-1]
    th = rng.normal(0.0, 0.6, circ.n_params())
    J, T, alpha = 0.7, 0.5, 0.18
    _, g, _, _, _, _ = B.hop_objective(circ, th, E, A, alpha, T, J)
    h = 1e-5
    TH = circ._shift_grid(th, h)
    vals = np.array([B.hop_objective(circ, TH[i], E, A, alpha, T, J)[0] for i in range(TH.shape[0])])
    g_fd = (vals[0::2] - vals[1::2]) / (2 * h)
    assert np.max(np.abs(g - g_fd)) < 1e-5
    assert _cos(g, g_fd) > 1 - 1e-8


def test_objective_is_not_a_function_of_p_for_J_positive():
    """Two states with IDENTICAL p and different amplitude signs differ in <psi|A|psi>."""
    rng = np.random.default_rng(5)
    dim = 32
    A = _rand_sym(dim, rng)
    psi = rng.normal(size=dim); psi /= np.linalg.norm(psi)
    psi2 = psi.copy(); psi2[::2] *= -1.0
    assert np.allclose(psi ** 2, psi2 ** 2)
    assert abs(B.hop_value(psi, A) - B.hop_value(psi2, A)) > 1e-3
    assert B.hop_abs_bound(psi, A) >= B.hop_value(psi, A) - 1e-12
    assert 0.0 <= B.sign_coherence(psi) <= 1.0
    assert B.sign_coherence(np.abs(psi)) == pytest.approx(1.0)


# ------------------------------------------------------------------------- the graph
def test_graph_symmetric_zero_diagonal_unit_spectral_norm_and_padding():
    rng = np.random.default_rng(1)
    W = rng.normal(0.0, 3.0, (60, 10, 3))
    D = B.pairwise_rmsd_matrix(W)
    assert np.allclose(D, D.T) and np.all(np.diag(D) == 0.0) and np.all(D >= 0)
    g = B.kernel_graph(D)
    A = g["A"]
    assert np.allclose(A, A.T)
    assert np.all(np.diag(A) == 0.0)
    assert np.linalg.eigvalsh(A)[-1] == pytest.approx(1.0, abs=1e-12)
    assert g["sigma"] == pytest.approx(float(np.median(D[np.triu_indices(60, 1)])))
    Ap = B.pad_graph(A, 64)
    assert Ap.shape == (64, 64)
    assert np.all(Ap[60:, :] == 0.0) and np.all(Ap[:, 60:] == 0.0)
    assert np.array_equal(Ap[:60, :60], A)
    # unit spectral norm bounds the hopping term by 1 for every unit vector
    for _ in range(5):
        v = rng.normal(size=64); v /= np.linalg.norm(v)
        assert B.hop_value(v, Ap) <= 1.0 + 1e-12


def test_perm_control_keeps_spectrum_and_degree_multiset():
    rng = np.random.default_rng(2)
    W = rng.normal(0.0, 3.0, (40, 9, 3))
    g = B.kernel_graph(B.pairwise_rmsd_matrix(W))
    perm = rng.permutation(40)
    gp = B.permuted_graph(g, perm)
    assert np.allclose(np.linalg.eigvalsh(gp["A_raw"]), np.linalg.eigvalsh(g["A_raw"]))
    assert np.allclose(np.sort(gp["degree"]), np.sort(g["degree"]))
    assert gp["lam_max"] == pytest.approx(g["lam_max"])
    assert not np.allclose(gp["degree"], g["degree"])          # the correspondence moved


def test_rand_control_matches_degree_sequence():
    rng = np.random.default_rng(4)
    W = rng.normal(0.0, 3.0, (80, 12, 3))
    g = B.kernel_graph(B.pairwise_rmsd_matrix(W))
    gr = B.random_degree_matched_graph(g, np.random.default_rng(9))
    A = gr["A_raw"]
    assert gr["sinkhorn_converged"]
    assert np.allclose(A, A.T) and np.all(np.diag(A) == 0.0) and np.all(A >= 0)
    assert np.max(np.abs(A.sum(1) - g["degree"])) / g["degree"].max() < 1e-8
    assert not np.allclose(A, g["A_raw"])                       # a different graph
    assert np.linalg.eigvalsh(gr["A"])[-1] == pytest.approx(1.0, abs=1e-12)


# ------------------------------------------------------------------- J = 0 reproduction
def test_J0_reproduces_run_cvar_vqe_bit_for_bit_synthetic():
    rng = np.random.default_rng(21)
    dim = 64
    E = np.concatenate([B.deployed_E(50), np.full(14, 12.0)])
    A = B.pad_graph(_rand_sym(50, rng), dim)
    for seed in (0, 1):
        p_ref, cv_ref, H_ref, _ = Q.run_cvar_vqe(E, 0.18, 0.5, n=6, layers=3, iters=20, seed=seed, lr=0.15)
        p, cv, Hn, hv, F, th, circ = B.run_hop_vqe(E, A, 0.0, alpha=0.18, T=0.5, n=6, layers=3,
                                                   iters=20, seed=seed, lr=0.15)
        assert np.array_equal(p, p_ref)
        assert cv == cv_ref and Hn == H_ref and hv == 0.0
        assert F == pytest.approx(cv_ref - 0.5 * H_ref)


@pytest.mark.skipif(not HAVE_CACHE, reason="S27 channel cache / S25 AMBER cache absent")
def test_J0_reproduces_arm_vqe_on_one_real_target():
    from s24 import d_harness as H
    from s27 import run_pool as RP
    from s22 import qcand_lib as QC
    cand, ch, _ = RP.channels_for("1A13")
    E = RP.zr(ch["DIS"])
    q = H.arm_vqe(cand, E, alpha=B.ALPHA, T=B.TEMP, layers=B.LAYERS, iters=B.ITERS, seed=0)
    enc = QC.Encoding(E)
    p, cv, Hn, hv, F, th, circ = B.run_hop_vqe(enc.E, np.zeros((enc.dim, enc.dim)), 0.0, seed=0, n=enc.n_qubits)
    r1 = B.readout_tail(enc, p, B.ALPHA)
    assert np.array_equal(r1["idx"], q["cands"])
    assert cv == q["cvar"] and Hn == q["entropy_nats"] and r1["m"] == q["m"]


# ------------------------------------------------------------------- the ground state
def test_ground_state_is_argmin_at_J0_and_eigenvector_at_J_positive():
    rng = np.random.default_rng(13)
    dim = 48
    E = B.deployed_E(dim)
    A = _rand_sym(dim, rng); A = A / np.linalg.eigvalsh(A)[-1]
    g0 = B.ground_state(E, A, 0.0)
    assert int(np.argmax(g0["p"])) == int(np.argmin(E))
    assert g0["p"][np.argmin(E)] == pytest.approx(1.0)
    g1 = B.ground_state(E, A, 1.0)
    Hm = np.diag(E) - 1.0 * A
    res = Hm @ g1["psi"] - g1["e0"] * g1["psi"]
    assert np.max(np.abs(res)) < 1e-10
    assert g1["p"].sum() == pytest.approx(1.0)
    assert g1["e0"] <= g0["e0"] + 1e-12                        # hopping can only lower E0
    assert 1.0 / np.sum(g1["p"] ** 2) > 1.0                    # delocalised at this J on a random A


# ---------------------------------------------------------------------- readouts and ties
def test_ptop_readout_breaks_ties_by_key_not_array_order():
    p = np.zeros(20); p[:3] = 0.3; p[3:] = 0.1 / 17            # 17 exact ties below the top 3
    key_a = np.random.default_rng(0).random(20)
    key_b = np.random.default_rng(1).random(20)
    ia = B.readout_ptop(p, 5, key_a)
    ib = B.readout_ptop(p, 5, key_b)
    assert set([0, 1, 2]) <= set(ia.tolist()) and set([0, 1, 2]) <= set(ib.tolist())
    assert not np.array_equal(ia, ib)                          # the key, not the order, decided


def test_weighted_readout_uniform_limit_equals_coordinate_average():
    from s12 import instrument as I
    rng = np.random.default_rng(17)
    W = rng.normal(0.0, 3.0, (30, 9, 3))
    D = B.pairwise_rmsd_matrix(W)
    C_w = B.readout_weighted(W, D, np.ones(30))
    C_u, _ = I.coordinate_average(W, P=D)
    assert np.allclose(C_w, C_u)


# ---------------------------------------------------------------------------- NaN-poison
def test_select_target_is_native_free_and_nan_poison_is_bit_identical(monkeypatch):
    """`select_target` never receives a native; `run_target` with NaN natives emits the same
    sets and coordinates and NaN RMSDs. A tiny synthetic pool keeps this under a few seconds."""
    from s24 import d_harness as H
    from s27 import run_pool as RP
    rng = np.random.default_rng(23)
    k, n = 40, 9
    W = rng.normal(0.0, 3.0, (k, n, 3))
    E = B.deployed_E(k)[rng.permutation(k)]
    key = RP.rng_for("SYN1", "tiekey").random(k)
    rows_a = B.select_target(W, E, key, "SYN1", j_grid=(0.0, 1.0), graphs=("REAL",), seeds=(0,))
    rows_b = B.select_target(W, E, key, "SYN1", j_grid=(0.0, 1.0), graphs=("REAL",), seeds=(0,))
    assert len(rows_a) == 4                                    # vqe J0, gs J0, vqe J1, gs J1
    for ra, rb in zip(rows_a, rows_b):
        for R in ("R1", "R3"):
            assert ra[R]["idx"] == rb[R]["idx"]
        for R in ("R1", "R2", "R3"):
            assert np.array_equal(np.asarray(ra[R]["C"]), np.asarray(rb[R]["C"]))
    # now through run_target with NaN natives: monkeypatch the channel loader
    nat = rng.normal(0.0, 3.0, (n, 3))
    cand_ok = H.Candidates(pdb="SYN1", n=n, seq="A" * n, fold=0, W=W, nat_ca=nat, oracle_rr=np.ones(k))
    cand_nan = H.Candidates(pdb="SYN1", n=n, seq="A" * n, fold=0, W=W, nat_ca=np.full((n, 3), np.nan),
                            oracle_rr=np.full(k, np.nan))
    calls = {"i": 0}

    def fake_channels(pdb):
        c = cand_ok if calls["i"] == 0 else cand_nan
        calls["i"] += 1
        return c, {"DIS": E}, {}

    monkeypatch.setattr(RP, "channels_for", fake_channels)
    monkeypatch.setattr(B, "J_GRID", (0.0, 1.0))
    monkeypatch.setattr(B, "GRAPHS", ("REAL",))
    monkeypatch.setattr(B, "SEEDS", (0,))
    ra = B.run_target("SYN1")
    rb = B.run_target("SYN1")
    assert len(ra) == len(rb) == 4
    for a, b in zip(ra, rb):
        for R in ("R1", "R2", "R3"):
            assert np.array_equal(np.asarray(a[R]["C"]), np.asarray(b[R]["C"]))
            assert np.isfinite(a[R]["rmsd"]) and np.isnan(b[R]["rmsd"])
        assert a["m"] == b["m"] and a["jac75"] == b["jac75"]
        assert np.isfinite(a["rmsd_dis75"]) and np.isnan(b["rmsd_dis75"])


def test_gate_still_passes_at_J_positive_by_construction():
    """The subset property is the tail-reading operator's, not the Hamiltonian's."""
    from s22 import qcand_lib as QC
    from s24 import d_harness as H
    rng = np.random.default_rng(29)
    k = 40
    E = B.deployed_E(k)[rng.permutation(k)]
    enc = QC.Encoding(E)
    A = B.pad_graph(_rand_sym(k, rng), enc.dim); A = A / np.linalg.eigvalsh(A)[-1]
    p, *_ = B.run_hop_vqe(enc.E, A, 3.0, n=enc.n_qubits, iters=15, seed=0)
    r1 = B.readout_tail(enc, p, B.ALPHA)
    gate = H.gate_set_equality(E, r1["idx"], r1["m"])
    assert gate["pass_"]
