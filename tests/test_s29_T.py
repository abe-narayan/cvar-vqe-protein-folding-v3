"""tests/test_s29_T.py -- lane T's derivations, pinned as executable checks.

Every test here is a property of the mathematics or of the code, not of the data: none reads a
native, an RMSD or a pool, and all of them run in under a second. They exist so that the claims in
`s29/THEORY.md` can be re-checked by another lane without reading the derivation.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from core import quantum as Q                      # noqa: E402
from s29 import s29_T_reach as R                   # noqa: E402
from s29 import s29_T_spectra as SP                # noqa: E402


# ------------------------------------------------------------------ section 3b: the variance law
def test_variance_law_matches_s28_hop_only_at_the_deployed_width():
    """THEORY 3b (3.4): Var[dF/dtheta] ~= r_stable(A)/D^2 for a unit-spectral-norm observable.

    Checked against a direct parameter-shift measurement on a synthetic near-rank-one graph, at
    the two widths S28 measured, to within a factor of 2 (the 2-design approximation's accuracy).
    """
    rng = np.random.default_rng(7)
    for n in (6, 8):
        D = 1 << n
        X = rng.normal(size=(D, 3))
        d = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
        sigma = float(np.median(d[np.triu_indices(D, 1)]))
        A = np.exp(-d ** 2 / (2 * sigma ** 2))
        np.fill_diagonal(A, 0.0)
        A = A / float(np.abs(np.linalg.eigvalsh(A)).max())
        r_stable = float((A ** 2).sum())            # ||A||_F^2 with ||A||_2 = 1
        meas = SP.__dict__  # keep the import used
        row = __import__("s27.s28_B_hop", fromlist=["x"]).measure_hop(
            n, 3, 0.18, 0.5, 1.0, 60, 11, 0.6, R.deployed_E(D), A, "hop_only")
        pred = r_stable / D ** 2
        assert 0.4 < row["var_g0"] / pred < 2.5, (n, row["var_g0"], pred, r_stable)
        assert meas is not None


def test_double_centering_removes_the_uniform_mode_exactly():
    """THEORY 3c: A_c = H A H has the uniform vector in its kernel to machine precision."""
    rng = np.random.default_rng(3)
    A = rng.normal(size=(64, 64))
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    Ac = SP.double_center(A)
    one = np.ones(64)
    assert np.abs(Ac @ one).max() < 1e-10
    w, v = np.linalg.eigh(Ac)
    top = v[:, int(np.argmax(np.abs(w)))]
    assert (top @ (one / np.sqrt(64))) ** 2 < 1e-20


def test_stable_rank_of_a_k_regular_graph_is_M_over_k():
    """THEORY 3b: the design rule's one exactly computable case (it predicts S28-B2's 46x)."""
    M, k = 500, 10
    A = np.zeros((M, M))
    for i in range(M):
        for j in range(1, k // 2 + 1):
            A[i, (i + j) % M] = 1.0
            A[(i + j) % M, i] = 1.0
    lam = float(np.abs(np.linalg.eigvalsh(A)).max())
    r_stable = float((A ** 2).sum()) / lam ** 2
    assert abs(r_stable - M / k) < 1e-6 * M / k


# ------------------------------------------------------- Q1: the flat set and the entropy optimum
def test_cvar_is_exactly_flat_above_the_var():
    """Q1.1: dCVaR/dp(x) = 0 for every state above the quantile (the 437-of-511 claim)."""
    D, alpha = 128, 0.18
    E = R.deployed_E(D)
    rng = np.random.default_rng(1)
    p = rng.random(D)
    p /= p.sum()
    val, q, dp = Q.cvar_exact(E, p, alpha)
    above = E > q
    assert above.sum() > 0
    assert np.abs(dp[above]).max() == 0.0
    assert (np.abs(dp[E < q]) > 0).all()


def test_simplex_optimum_reduces_to_the_gibbs_state_at_alpha_one():
    """Q1.2: the Rockafellar-Uryasev solution must return -T log Z and the Boltzmann law."""
    E = R.deployed_E(128)
    T = 0.3
    out = R.simplex_optimum(E, 1.0, T)
    gibbs = float(-T * np.log(np.exp(-E / T).sum()))
    assert abs(out["F"] - gibbs) < 1e-6
    p = np.exp(-E / T)
    p /= p.sum()
    H = float(-(p * np.log(p)).sum())
    assert abs(out["entropy_nats"] - H) < 1e-6


def test_simplex_optimum_is_attained_and_beats_the_uniform_state():
    """Q1.2: F(p*) evaluated directly must equal the RU value, and must be below F(uniform)."""
    D, alpha, T = 512, 0.18, 0.5
    E = R.deployed_E(D)
    out = R.simplex_optimum(E, alpha, T)
    c = np.maximum(out["t_star"] - E, 0.0) / alpha
    w = np.exp((c - c.max()) / T)
    ps = w / w.sum()
    cv, _, _ = Q.cvar_exact(E, ps, alpha)
    F_direct = cv - T * float(-(ps * np.log(ps)).sum())
    # the outer maximisation is on a 20,001-point grid over a range of about 7, and F is flat to
    # second order at its optimum, so the RU value and the direct evaluation agree to the grid.
    assert abs(F_direct - out["F"]) < 1e-4
    pu = np.full(D, 1.0 / D)
    cvu, _, _ = Q.cvar_exact(E, pu, alpha)
    F_uniform = cvu - T * float(-(pu * np.log(pu)).sum())
    assert out["F"] < F_uniform


# ------------------------------------------------------ section 4.5: the set-equality counterexample
def test_aggregate_objective_breaks_the_prefix_rule():
    """THEORY 4.5: the optimal 2-subset is the two WORST states by every per-state criterion."""
    W = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 0.1]])
    t = np.zeros(2)
    per_state = np.linalg.norm(W - t, axis=1)
    order = np.argsort(per_state)                       # 3 is best individually
    assert order[0] == 2
    f_prefix = min(np.linalg.norm(W[[order[0], order[1]]].mean(0) - t),
                   np.linalg.norm(W[[order[0], order[2]]].mean(0) - t))
    f_best = np.linalg.norm(W[[0, 1]].mean(0) - t)
    assert f_best < f_prefix - 0.4
    assert f_best == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------- section 1.2: the contraction identity
def test_contraction_is_an_exact_variance_identity():
    """THEORY 1.2 (1.2)/(1.3): d(mean)^2 = mean d_k^2 - spread^2, exactly."""
    rng = np.random.default_rng(5)
    K, N = 40, 12
    W = rng.normal(size=(K, N, 3))
    C = W.mean(0)
    i, j = 0, 5
    v = W[:, i, :] - W[:, j, :]
    d_k = np.linalg.norm(v, axis=1)
    s2 = float((np.linalg.norm(v - v.mean(0), axis=1) ** 2).mean())
    lhs = float(np.linalg.norm(C[i] - C[j])) ** 2
    assert abs(lhs - ((d_k ** 2).mean() - s2)) < 1e-10
    assert lhs <= (d_k.mean()) ** 2 + 1e-12          # Jensen, the weaker corollary


# --------------------------------------------------------- section 5.3: the prefix cut-rank asymmetry
def test_prefix_sets_are_low_rank_and_random_sets_are_not():
    """THEORY 5.3: the indicator of a prefix has cut rank <= 2; a random 75-subset saturates."""
    n = 9
    D = 1 << n
    v = np.zeros(D)
    v[:75] = 1.0
    prefix_rank = max(np.linalg.matrix_rank(v.reshape(1 << k, 1 << (n - k))) for k in range(1, n))
    rng = np.random.default_rng(0)
    v2 = np.zeros(D)
    v2[rng.choice(D, 75, replace=False)] = 1.0
    rand_rank = max(np.linalg.matrix_rank(v2.reshape(1 << k, 1 << (n - k))) for k in range(1, n))
    assert prefix_rank <= 2
    assert rand_rank >= 12


# ------------------------------------------------------------------ Q2: the local mixer's stable rank
def test_transverse_field_stable_rank_is_D_over_n():
    """Q2 (Q2.2): r_stable(sum_q X_q) = D/n, the one non-commuting operator that scales."""
    for n in (4, 6, 8):
        D = 1 << n
        M = np.zeros((D, D))
        for q in range(n):
            bit = 1 << (n - 1 - q)
            idx = np.arange(D)
            M[idx, idx ^ bit] += 1.0
        fro2 = float((M ** 2).sum())
        lam = float(np.abs(np.linalg.eigvalsh(M)).max())
        assert abs(fro2 - n * D) < 1e-8
        assert abs(lam - n) < 1e-8
        assert abs(fro2 / lam ** 2 - D / n) < 1e-8


# ------------------------------------------------------------ section 8: the displacement bound
def test_displacement_identity_and_the_cosine_price_table():
    """THEORY 8.1 (8.2): RMSD(c + rho*|e|*u_hat) = RMSD_prod * sqrt(1 - rho^2), and the price
    table (8.3) inverts it. Checked as exact algebra in the aligned frame."""
    rng = np.random.default_rng(11)
    N = 12
    t = rng.normal(size=(N, 3))
    c = t + rng.normal(size=(N, 3))
    e = (t - c).ravel()
    ne = float(np.linalg.norm(e))
    for rho in (0.04, 0.14, 0.37, 0.628):
        # a displacement with exactly this cosine, at the optimal step size s = rho
        g = rng.normal(size=e.shape)
        g -= (g @ e) / ne ** 2 * e
        g /= np.linalg.norm(g)
        u_hat = rho * e / ne + np.sqrt(1 - rho ** 2) * g
        u = rho * ne * u_hat
        got = float(np.linalg.norm(c.ravel() + u - t.ravel()))
        assert abs(got - ne * np.sqrt(1 - rho ** 2)) < 1e-9
    prod = 3.2126
    assert abs(prod * np.sqrt(1 - 0.628 ** 2) - 2.50) < 0.005      # the 2.5 A price
    assert abs(prod * np.sqrt(1 - 0.140 ** 2) - 3.181) < 0.002     # the random-field price
    assert abs(prod * np.sqrt(1 - 0.374 ** 2) - 2.98) < 0.005      # PC1 with a perfect sign


def test_sign_accuracy_enters_as_two_q_minus_one():
    """THEORY 8.2 (8.4): a field with per-target magnitude |rho| and sign accuracy q enters at
    |rho|(2q-1); at PC1's measured 52% that is 4% of its magnitude."""
    rho, q = 0.37, 0.52
    assert abs(rho * (2 * q - 1) - 0.0148) < 1e-4
    gain = 3.2126 * (1 - np.sqrt(1 - (rho * (2 * q - 1)) ** 2))
    assert gain < 0.001                                            # indistinguishable from zero
