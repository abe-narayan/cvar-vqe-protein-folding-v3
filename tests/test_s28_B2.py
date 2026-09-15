"""Tests for `s27/s28_B2_knn.py` (S28 lane B2: the kNN spread-spectrum hopping graph).

Pins: the kNN graph is symmetric with zero diagonal, unit spectral norm, every degree >= k
(symmetrisation can only add edges), zero padding rows; a degree-normalised connected graph has
lambda_1 = 1 before the rescale; the spectral report puts the kNN numbers beside the Gaussian
graph's and, on a pool with cluster structure, the kNN spectrum is spread (lambda_2 / lambda_1
well above the Gaussian graph's); the hopping gradient on the kNN graph matches finite
differences; the graph switch in `run_main` is restored after use.

Run:  python -m pytest tests/test_s28_B2.py -q -p no:cacheprovider
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from core import quantum as Q                      # noqa: E402
from s27 import s28_B_hop as B                     # noqa: E402
from s27 import s28_B2_knn as B2                   # noqa: E402


def _pool(seed=0, k=60, n=9, clusters=4):
    rng = np.random.default_rng(seed)
    centres = rng.normal(0.0, 6.0, (clusters, n, 3))
    W = np.concatenate([c[None] + rng.normal(0.0, 1.0, (k // clusters, n, 3)) for c in centres], 0)
    return W


def test_knn_graph_symmetric_zero_diag_unit_norm_degree_at_least_k():
    W = _pool()
    D = B.pairwise_rmsd_matrix(W)
    for k in (3, 5, 10):
        g = B2.knn_graph(D, k)
        A = g["A"]
        assert np.allclose(A, A.T)
        assert np.all(np.diag(A) == 0.0)
        assert np.linalg.eigvalsh(A)[-1] == pytest.approx(1.0, abs=1e-12)
        assert g["knn_degree_min"] >= k
        assert np.all(g["A_bin"] == g["A_bin"].T)
        assert set(np.unique(g["A_bin"]).tolist()) <= {0.0, 1.0}
        # each row lists at least its own k nearest
        Dm = D.copy(); np.fill_diagonal(Dm, np.inf)
        for i in range(len(D)):
            nn = np.argsort(Dm[i], kind="stable")[:k]
            assert np.all(g["A_bin"][i, nn] == 1.0)
    Ap = B.pad_graph(B2.knn_graph(D, 5)["A"], 64)
    assert np.all(Ap[60:, :] == 0.0) and np.all(Ap[:, 60:] == 0.0)


def test_degree_normalised_connected_graph_has_lambda1_one_before_rescale():
    W = _pool(seed=3)
    D = B.pairwise_rmsd_matrix(W)
    g = B2.knn_graph(D, 10)
    if g["n_components"] == 1:
        assert g["lam_max"] == pytest.approx(1.0, abs=1e-10)
    else:
        assert g["lam_max"] <= 1.0 + 1e-10


def test_knn_spectrum_is_spread_relative_to_gaussian_on_a_clustered_pool():
    W = _pool(seed=5, k=80, clusters=4)
    D = B.pairwise_rmsd_matrix(W)
    rep = B2.spectral_report(D, 5)
    assert rep["knn"]["lam2_over_lam1"] > rep["gauss"]["lam2_over_lam1"]
    assert 0.0 < rep["knn"]["perron_pr_over_dim"] <= 1.0
    assert rep["knn"]["density"] < 0.5


def test_hop_gradient_on_knn_graph_matches_finite_differences():
    rng = np.random.default_rng(11)
    circ = Q.StatevectorCircuit(6, 3)
    W = _pool(seed=7, k=50)
    A = B.pad_graph(B2.knn_graph(B.pairwise_rmsd_matrix(W), 5)["A"], circ.dim)
    th = rng.normal(0.0, 0.6, circ.n_params())
    assert np.max(np.abs(B.grad_hop_paramshift(circ, th, A) - B.grad_hop_fd(circ, th, A))) < 1e-6


def test_k_out_of_range_raises():
    D = B.pairwise_rmsd_matrix(_pool(k=20))
    with pytest.raises(ValueError):
        B2.knn_graph(D, 0)
    with pytest.raises(ValueError):
        B2.knn_graph(D, 20)


def test_run_main_restores_the_gaussian_graph_switch(monkeypatch):
    """The endpoint runner swaps `s28_B_hop.kernel_graph` for the kNN builder and must put it
    back even when the run raises."""
    original = B.kernel_graph

    def boom(pdb):
        raise RuntimeError("stop")

    monkeypatch.setattr(B, "run_target", boom)
    from s25 import phys_lib as P
    monkeypatch.setattr(P, "targets", lambda: ["XXXX"])
    with pytest.raises(RuntimeError):
        B2.run_main(10, limit=1)
    assert B.kernel_graph is original
