"""SPRINT 14 / VQE -- correctness tests for the CVaR estimator and its gradient.

THE DERIVATION, WRITTEN DOWN
============================
Everything below tests one derivation, so it is stated once here and each test names the
line it checks.

Let ``p_theta`` be the distribution the circuit prepares over basis states ``x``, and
``E(x)`` the objective.  We want the mean of the LOWER ``alpha`` tail of the energy, which
is the Rockafellar-Uryasev variational form with the sign flipped for a minimisation:

    (1)    CVaR_alpha(theta) = max_t [ t - (1/alpha) E_{x~p}[ (t - E(x))_+ ] ]

The inner function ``G(t) = t - (1/alpha) E[(t - E)_+]`` is concave with

    (2)    dG/dt = 1 - (1/alpha) P(E < t),

so it is maximised where ``P(E < t) = alpha``, i.e. at the alpha-quantile ``q``.  Substituting
back, and being careful that a DISCRETE distribution generally has no state at which the
cumulative mass equals alpha exactly, the boundary atom is split fractionally:

    (3)    CVaR_alpha = (1/alpha) [ sum_{E(x) < q} p(x) E(x) + (alpha - P(E < q)) q ] .

That is exactly what `core.quantum.cvar_exact` computes, and equation (3) -- not "the mean
of the lowest ceil(alpha N) samples" -- is the object being differentiated.

Differentiating (1) at the optimal ``t``, the envelope theorem kills the ``t`` dependence
(``dG/dt = 0`` there), leaving only the explicit ``theta`` dependence:

    (4)    grad CVaR = -(1/alpha) sum_x grad p(x) (q - E(x))_+
                     = E_p[ w(X) grad log p(X) ] ,    w(x) = (E(x) - q)/alpha if E(x) < q
                                                             0                otherwise.

``w`` is exactly ``dCVaR/dp``, which is why `cvar_exact` returns it as ``dp``.

THE BASELINE, AND THE EXACT FORM OF THE RECORDED DEFECT
-------------------------------------------------------
For any CONSTANT ``c``, ``E_p[c grad log p] = c grad sum_x p(x) = c grad 1 = 0``, so
``w -> w - c`` changes the variance and not the mean.  That is the ``baseline="const"``
form and it is unbiased.

The recorded defect subtracts a constant from the TAIL ENTRIES ONLY, i.e. it uses
``b(x) = c * 1[x in tail]``, which is a FUNCTION OF x and not a constant.  Its expectation
is not zero, and the bias has a closed form:

    (5)    E[g_tail] - grad CVaR = - c * grad_theta P(E(X) < q) ,

with ``c`` the p-weighted mean of ``w`` over the tail.  So the defect does not merely add
noise: it adds ``c`` times the gradient of the TAIL PROBABILITY, an entirely different
direction that pushes the optimiser to change how much mass sits in the tail rather than
how good the tail is.  `t_tail_baseline_bias_has_the_closed_form_of_equation_5` checks (5)
to machine precision, which is stronger than measuring a cosine.

THE ZERO-GRADIENT REGIME
------------------------
From (3): if the single lowest-energy state already carries mass ``>= alpha`` then the sum
over ``E < q`` is empty, ``q = min E``, and ``CVaR_alpha = min E`` identically.  The
gradient is then EXACTLY zero -- not small, not a barren plateau, but an exact flat region
of measure-zero-free extent.  It is also the global minimum of CVaR_alpha, so it is
convergence and not failure; but it means that at small alpha the CVaR objective is a pure
argmin-finder with no pressure on anything else, which is the mechanism behind the recorded
trap that "at alpha <= 0.25 two AMBER variants became literally the same objective".
`t_cvar_gradient_is_exactly_zero_iff_pmin_ge_alpha` states and checks the iff.

    python -m pytest s14/test_vqe.py -q
"""
from __future__ import annotations

import os
import sys
from decimal import Decimal, getcontext

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from core import quantum as Q            # noqa: E402
import qansatz                           # noqa: E402


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0.0 or nb == 0.0:
        return float("nan")
    return float(a @ b / (na * nb))


def _brute_cvar(E, p, alpha):
    """Equation (3) evaluated by walking the mass in 60-digit decimal. Independent code."""
    getcontext().prec = 60
    E = np.asarray(E, float)
    p = np.asarray(p, float)
    p = p / p.sum()
    a = Decimal(repr(float(alpha)))
    acc = Decimal(0)
    tot = Decimal(0)
    for i in np.argsort(E, kind="stable"):
        take = min(Decimal(repr(float(p[i]))), a - acc)
        if take <= 0:
            break
        tot += take * Decimal(repr(float(E[i])))
        acc += take
        if acc >= a:
            break
    return float(tot / acc)


# ===================================================================== the value
def test_cvar_exact_is_equation_3_including_the_partial_bucket():
    """400 random discrete distributions against a 60-digit independent reference."""
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(400):
        m = int(rng.integers(2, 40))
        E = rng.normal(0, 1, m)
        p = rng.random(m) ** 2 + 1e-6
        p /= p.sum()
        a = float(rng.choice([1.0, 0.5, 0.37, 0.25, 0.1, 0.05, 0.01, 0.003]))
        worst = max(worst, abs(Q.cvar_exact(E, p, a)[0] - _brute_cvar(E, p, a)))
    assert worst < 1e-13, worst


def test_two_independent_implementations_agree():
    """`cvar_exact` and `cvar_from_probs` are separate code paths; they must coincide."""
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(300):
        m = int(rng.integers(2, 60))
        E = rng.normal(0, 1, m)
        p = rng.random(m) + 1e-9
        p /= p.sum()
        a = float(rng.uniform(0.005, 1.0))
        worst = max(worst, abs(Q.cvar_exact(E, p, a)[0] - Q.cvar_from_probs(E, p, a)[0]))
    assert worst < 1e-12, worst


def test_ties_do_not_leak_index_order_into_the_value():
    """With degenerate energies the answer must not depend on how the states are ordered."""
    rng = np.random.default_rng(2)
    E = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 2.0])
    p = np.array([0.3, 0.1, 0.05, 0.2, 0.25, 0.1])
    for a in (0.05, 0.2, 0.4, 0.45, 0.7, 1.0):
        base = Q.cvar_exact(E, p, a)[0]
        for _ in range(60):
            perm = rng.permutation(len(E))
            assert abs(Q.cvar_exact(E[perm], p[perm], a)[0] - base) < 1e-12
        assert abs(base - _brute_cvar(E, p, a)) < 1e-13


def test_massive_degeneracy_is_the_fractional_answer_not_the_bucket_answer():
    """24 states on 3 distinct energies: the boundary bucket must split, not round."""
    E = np.repeat(np.array([0.0, 1.0, 2.0]), 8).astype(float)
    p = np.full(24, 1 / 24)
    # alpha = 0.4 takes all of the 8/24 at E=0 plus (0.4 - 1/3) of the E=1 block
    v = Q.cvar_exact(E, p, 0.4)[0]
    expect = ((1 / 3) * 0.0 + (0.4 - 1 / 3) * 1.0) / 0.4
    assert abs(v - expect) < 1e-13, (v, expect)


def test_alpha_below_the_largest_atom_returns_the_minimum_over_the_support():
    """The ``alpha * 2^n < 1`` edge case -- and a subtlety worth stating explicitly.

    With mass on only two states the value is the minimum over the SUPPORT, for every
    alpha at or below the modal mass, including ``alpha < 2^-n``.

    But ``dCVaR/dp`` is NOT zero here, and that is correct rather than a defect: a state
    with ``p(x) = 0`` and ``E(x) < q`` has a perfectly good one-sided derivative
    ``(E(x) - q)/alpha < 0``, meaning "moving mass onto this state would lower the CVaR".
    This is precisely the term that lets a variational optimiser DISCOVER good states it is
    not yet sampling, so the zero-gradient regime is governed by the argmin of ``E`` over
    the whole register and not by the argmin over the current support.
    """
    rng = np.random.default_rng(3)
    n = 6
    E = rng.normal(0, 1, 1 << n)
    p = np.zeros(1 << n)
    p[3], p[7] = 0.9, 0.1
    for a in (0.5, 0.1, 0.05, 1.0 / (1 << n), 1e-6):
        v, q, dp = Q.cvar_exact(E, p, a)
        assert abs(v - min(E[3], E[7])) < 1e-13
        # every strictly-better state that is currently unoccupied pulls, and only those
        assert np.array_equal(dp != 0.0, E < q)
        assert np.all(dp[E < q] < 0.0)
    # with STRICTLY POSITIVE p the zero-gradient regime is recovered exactly
    q2 = np.full(1 << n, 1e-9)
    q2[int(np.argmin(E))] = 1.0
    q2 /= q2.sum()
    assert not np.any(Q.cvar_exact(E, q2, 0.5)[2] != 0.0)


# =============================================================== the zero regime
def test_cvar_gradient_is_exactly_zero_iff_pmin_ge_alpha():
    """The iff from the docstring: 3,000 random cases, zero counterexamples permitted."""
    bad = []
    for t in range(3000):
        rng = np.random.default_rng(t)
        m = int(rng.integers(2, 30))
        E = rng.normal(0, 1, m)
        p = rng.random(m) ** 3 + 1e-9
        p /= p.sum()
        a = float(rng.choice([0.5, 0.25, 0.1, 0.05, 0.01]))
        _, _, dp = Q.cvar_exact(E, p, a)
        if (p[int(np.argmin(E))] >= a) != (not np.any(dp != 0)):
            bad.append(t)
    assert not bad, bad[:10]


def test_the_zero_regime_is_the_global_minimum_of_cvar_not_a_plateau():
    """CVaR_alpha >= min E always, with equality exactly in the zero-gradient regime."""
    rng = np.random.default_rng(4)
    for _ in range(400):
        m = int(rng.integers(2, 40))
        E = rng.normal(0, 1, m)
        p = rng.random(m) ** 3 + 1e-9
        p /= p.sum()
        a = float(rng.choice([0.5, 0.2, 0.05, 0.01]))
        v = Q.cvar_exact(E, p, a)[0]
        assert v >= E.min() - 1e-12
        flat = not np.any(Q.cvar_exact(E, p, a)[2] != 0)
        assert flat == (abs(v - E.min()) < 1e-12)


# ================================================================= the gradient
GRID = [(n, L, a) for n in (5, 7, 9) for L in (1, 2, 3, 5)
        for a in (1.0, 0.5, 0.25, 0.1, 0.05, 0.01)]


@pytest.mark.parametrize("n,layers,alpha", GRID)
def test_paramshift_equals_finite_differences(n, layers, alpha):
    """Equation (4) against wholly independent machinery, over alpha x n x depth x seed.

    Finite differences on the exact CVaR VALUE use no projector algebra and no envelope
    theorem, so agreement is a real check on the analytic gradient rather than a
    restatement of it.  Cells whose gradient is exactly zero (the regime above) are
    checked for that instead, since a cosine is undefined there.
    """
    for seed in range(4):
        rng = np.random.default_rng(1000 * seed + 13 * n + layers)
        circ = Q.StatevectorCircuit(n, layers)
        th = rng.normal(0, 0.8, circ.n_params())
        E = rng.normal(0, 1, 1 << n)
        gp = Q.grad_cvar_paramshift(circ, th, E, alpha)
        gf = Q.grad_cvar_fd(circ, th, E, alpha, h=1e-5)
        if np.linalg.norm(gp) == 0.0:
            assert np.linalg.norm(gf) == 0.0
            assert Q.cvar_exact(E, circ.probs(th), alpha)[0] == pytest.approx(E.min())
            continue
        assert _cos(gp, gf) > 1 - 1e-12, _cos(gp, gf)
        assert np.linalg.norm(gp - gf) / np.linalg.norm(gf) < 1e-6


def test_paramshift_equals_exact_differentiation_of_cvar_exact():
    """Chain rule check: grad CVaR = (dCVaR/dp) . (dp/dtheta), both computed separately."""
    rng = np.random.default_rng(5)
    for n, layers, alpha in [(6, 2, 0.3), (8, 3, 0.1), (7, 1, 0.5)]:
        circ = Q.StatevectorCircuit(n, layers)
        th = rng.normal(0, 0.9, circ.n_params())
        E = rng.normal(0, 1, 1 << n)
        dp = Q.cvar_exact(E, circ.probs(th), alpha)[2]
        # dp/dtheta by central differences on the probabilities alone
        h = 1e-6
        J = np.empty((circ.n_params(), 1 << n))
        for kk in range(circ.n_params()):
            tp = th.copy(); tp[kk] += h
            tm = th.copy(); tm[kk] -= h
            J[kk] = (circ.probs(tp) - circ.probs(tm)) / (2 * h)
        assert _cos(J @ dp, Q.grad_cvar_paramshift(circ, th, E, alpha)) > 1 - 1e-9


# ================================================================== the baseline
def test_const_baseline_is_unbiased_and_tail_baseline_is_not():
    """Exact EXPECTATIONS of both estimators -- zero sampling noise, so any gap is bias."""
    rng = np.random.default_rng(20260905)
    for alpha in (0.5, 0.25, 0.1, 0.05):
        n, layers = 7, 3
        th = rng.normal(0, 0.8, n * layers)
        E = rng.normal(0, 1.0, 1 << n)
        circ = Q.StatevectorCircuit(n, layers)
        p = circ.probs(th)
        g_true = Q.grad_cvar_paramshift(circ, th, E, alpha)
        q = Q.cvar_exact(E, p, alpha)[1]
        PR = circ.probs_batch(circ._shift_grid(th, np.pi / 2))
        G = ((PR[0::2] - PR[1::2]) / 2.0 / np.maximum(p, 1e-15)[None, :]).T
        w = np.where(E < q, (E - q) / alpha, 0.0)
        g_const = (p[:, None] * (w - (p * w).sum())[:, None] * G).sum(0)
        m = E < q
        wt = np.zeros_like(w)
        wt[m] = w[m] - (p[m] * w[m]).sum() / p[m].sum()
        g_tail = (p[:, None] * wt[:, None] * G).sum(0)
        assert _cos(g_const, g_true) > 1 - 1e-9
        assert abs(np.linalg.norm(g_const) / np.linalg.norm(g_true) - 1) < 1e-9
        assert _cos(g_tail, g_true) < 0.99


def test_tail_baseline_bias_has_the_closed_form_of_equation_5():
    """The defect is not noise: its bias is exactly ``-c * grad P(E < q)``.

    Pinning the MECHANISM is stronger than pinning a cosine, because a cosine varies with
    alpha and instance (measured here from +0.32 to +0.96) while the formula does not.
    """
    rng = np.random.default_rng(77)
    for alpha in (0.5, 0.25, 0.1, 0.05):
        n, layers = 7, 3
        th = rng.normal(0, 0.8, n * layers)
        E = rng.normal(0, 1.0, 1 << n)
        circ = Q.StatevectorCircuit(n, layers)
        p = circ.probs(th)
        q = Q.cvar_exact(E, p, alpha)[1]
        PR = circ.probs_batch(circ._shift_grid(th, np.pi / 2))
        dP = (PR[0::2] - PR[1::2]) / 2.0                    # (P, dim) = dp(x)/dtheta_k
        G = (dP / np.maximum(p, 1e-15)[None, :]).T
        w = np.where(E < q, (E - q) / alpha, 0.0)
        m = E < q
        c = (p[m] * w[m]).sum() / p[m].sum()
        wt = np.zeros_like(w)
        wt[m] = w[m] - c
        g_tail = (p[:, None] * wt[:, None] * G).sum(0)
        g_true = Q.grad_cvar_paramshift(circ, th, E, alpha)
        predicted_bias = -c * dP[:, m].sum(1)               # -c * grad P(E < q)
        assert np.abs(g_tail - g_true - predicted_bias).max() < 1e-9, alpha


def test_qansatz_cvar_gradient_is_the_defect_bit_for_bit():
    """The shipped `qansatz` estimator IS `baseline="tail"`; the defect is reproduced."""
    rng = np.random.default_rng(11)
    n = 10
    an_new = Q.OneLayerAnsatz(n, ring=True)
    an_old = qansatz.OneLayerAnsatz(n, ring=True)
    th = rng.uniform(0.3, 2.8, n)
    for alpha in (0.5, 0.25, 0.1, 0.05):
        b = an_new.sample(th, 4096, np.random.default_rng(3))
        e = rng.standard_normal(4096)
        g_core = Q.cvar_gradient(an_new, th, b, e, alpha, baseline="tail")[0]
        g_old = qansatz.cvar_gradient(an_old, th, b, e, alpha)[0]
        assert np.abs(g_core - g_old).max() == 0.0, alpha


def test_sampled_const_baseline_converges_to_the_exact_gradient():
    """The deployable estimator, as a sampled quantity, must have the right mean."""
    n = 10
    rng = np.random.default_rng(12)
    an = Q.OneLayerAnsatz(n, ring=True)
    th = rng.uniform(0.3, 2.8, n)
    bits_all = Q.all_bitstrings(n)
    E = rng.normal(0, 1, 1 << n)
    for alpha in (0.5, 0.25, 0.1):
        g_ex = Q.cvar_gradient_exact(an, th, bits_all, E, alpha)[0]
        acc = []
        for s in range(8):
            b = an.sample(th, 40000, np.random.default_rng(s))
            idx = b.astype(np.int64) @ (1 << np.arange(n - 1, -1, -1))
            acc.append(Q.cvar_gradient(an, th, b, E[idx], alpha, baseline="const")[0])
        g = np.mean(acc, axis=0)
        assert _cos(g, g_ex) > 0.999, (alpha, _cos(g, g_ex))
        assert abs(np.linalg.norm(g) / np.linalg.norm(g_ex) - 1) < 0.02


# =================================================== the sampled value estimator
def test_ceil_cvar_estimator_is_biased_when_alpha_times_N_is_not_an_integer():
    """`cvar_from_samples` averages the lowest ceil(alpha N); equation (3) splits the atom.

    They coincide exactly when ``alpha * N`` is an integer and differ otherwise, with a
    bias that is UPWARD (the extra state is worse than the tail it joins) and decays like
    1/N.  This is not a defect to fix -- it is the estimator the SPSA driver evaluates --
    but it must be known, because at small alpha and modest shot counts it is large.
    """
    rng = np.random.default_rng(13)
    for N, a, lo in [(13, 0.1, 0.05), (10, 0.25, 0.03), (8, 0.3, 0.04)]:
        d = [Q.cvar_from_samples(e, a) - Q.cvar_from_probs(e, np.full(N, 1 / N), a)[0]
             for e in rng.normal(0, 1, (3000, N))]
        assert np.mean(d) > lo, (N, a, np.mean(d))
    # integer alpha*N: exact agreement
    for N, a in [(100, 0.25), (1000, 0.1), (64, 0.5)]:
        e = rng.normal(0, 1, N)
        assert abs(Q.cvar_from_samples(e, a)
                   - Q.cvar_from_probs(e, np.full(N, 1 / N), a)[0]) < 1e-12
    # and the bias shrinks like 1/N at fixed alpha
    b = []
    for N in (16, 64, 256, 1024):
        d = [Q.cvar_from_samples(e, 0.1) - Q.cvar_from_probs(e, np.full(N, 1 / N), 0.1)[0]
             for e in rng.normal(0, 1, (2000, N))]
        b.append(abs(np.mean(d)))
    assert b[0] > b[1] > b[2] > b[3]


def test_tail_mask_matches_a_stable_sort_at_every_size():
    """`tail_indices` switches from argsort to argpartition at 448; the mask must not."""
    rng = np.random.default_rng(14)
    for N in (10, 447, 448, 449, 1000, 5000):
        for a in (1.0, 0.5, 0.25, 0.1, 0.013):
            e = rng.integers(0, 5, N).astype(float)      # heavy ties on purpose
            k, q, tail = Q.tail_indices(e, a)
            order = np.argsort(e, kind="stable")
            ref = np.zeros(N, bool)
            ref[order[:k]] = True
            assert np.array_equal(tail, ref), (N, a)
            assert q == e[order[k - 1]]


# ================================================= alpha collapses the objective
def test_small_alpha_collapses_distinct_objectives_onto_each_other():
    """The recorded trap, given its mechanism.

    At small alpha the CVaR of a concentrated distribution reads only the very bottom of
    the objective, so two objectives that agree on their lowest states are LITERALLY the
    same CVaR function even though they differ everywhere else.  Constructed here as an
    identity rather than measured as a coincidence: two objectives sharing their argmin
    give bit-identical CVaR and bit-identical gradient for every alpha at or below the
    modal mass.
    """
    rng = np.random.default_rng(15)
    n, layers = 7, 3
    circ = Q.StatevectorCircuit(n, layers)
    th = rng.normal(0, 0.4, circ.n_params())              # deliberately concentrated
    p = circ.probs(th)
    mode = int(np.argmax(p))
    E1 = rng.normal(0, 1, 1 << n)
    E1[mode] = E1.min() - 1.0                             # the modal state IS the argmin
    E2 = E1 + rng.uniform(1.0, 5.0, len(E1))              # differ everywhere ...
    E2[mode] = E1[mode]                                   # ... except at the argmin
    alpha = 0.9 * p.max()                                 # alpha <= p(argmin): collapse
    v1, q1, d1 = Q.cvar_exact(E1, p, alpha)
    v2, q2, d2 = Q.cvar_exact(E2, p, alpha)
    assert v1 == v2 == E1[mode] and np.array_equal(d1, d2)
    assert np.array_equal(Q.grad_cvar_paramshift(circ, th, E1, alpha),
                          Q.grad_cvar_paramshift(circ, th, E2, alpha))
    # ... and at alpha = 1 they are emphatically NOT the same objective
    assert Q.cvar_exact(E1, p, 1.0)[0] != Q.cvar_exact(E2, p, 1.0)[0]
    assert not np.array_equal(Q.grad_cvar_paramshift(circ, th, E1, 1.0),
                              Q.grad_cvar_paramshift(circ, th, E2, 1.0))


def test_alpha_schedule_is_monotone_and_hits_its_endpoints():
    assert Q.alpha_schedule(0.0, 0.5, 0.05) == pytest.approx(0.5)
    assert Q.alpha_schedule(1.0, 0.5, 0.05) == pytest.approx(0.05)
    xs = [Q.alpha_schedule(t / 20, 0.5, 0.05) for t in range(21)]
    assert all(xs[i] > xs[i + 1] for i in range(20))
    assert Q.alpha_schedule(-3.0, 0.5, 0.05) == pytest.approx(0.5)
    assert Q.alpha_schedule(9.0, 0.5, 0.05) == pytest.approx(0.05)


def test_invalid_alpha_is_rejected_everywhere():
    e = np.arange(10.0)
    p = np.full(10, 0.1)
    for a in (0.0, -0.1, 1.5, np.nan):
        for fn, args in ((Q.cvar, (e, a)), (Q.cvar_from_samples, (e, a)),
                         (Q.cvar_from_probs, (e, p, a)), (Q.tail_indices, (e, a))):
            with pytest.raises(ValueError):
                fn(*args)


# ============================================ the budget-honest driver's accounting
def test_the_evaluation_counter_is_a_hard_cap():
    from s14 import vqe_lib as V
    E = np.random.default_rng(0).standard_normal(4096)
    for budget in (1, 7, 100, 4095):
        c = V.Counter(E, budget)
        for _ in range(50):
            c(np.random.default_rng(0).integers(0, 4096, 137))
        assert c.used == budget, (budget, c.used)
        assert c.all_seen().size == budget


def test_every_classical_arm_respects_the_budget_exactly():
    from s14 import vqe_lib as V
    E = np.random.default_rng(1).standard_normal(4 ** 7)
    for name, fn in V.CLASSICAL.items():
        for budget in (200, 1500):
            c = fn(E, 7, 4, budget, np.random.default_rng(0))
            assert c.used == budget, (name, budget, c.used)
            assert E[c.best_i] == c.best_e


def test_vqe_run_spends_exactly_its_budget_and_reports_its_initialisation():
    from s14 import vqe_run as R
    E = np.random.default_rng(2).standard_normal(1 << 12)
    rmsd = np.random.default_rng(3).uniform(1, 6, 1 << 12)
    r = R.run(E, 12, 0.25, budget=2048, shots=256, ansatz="mps2f", seed=0, rmsd=rmsd)
    assert r["evals"] == 2048 and r["iters"] == 8
    assert np.isfinite(r["init_rmsd_mean"]) and np.isfinite(r["rmsd_returned"])
    assert 0.0 <= r["entropy_bits"] <= 12.0 + 1e-9
    assert r["baseline"] == "const"
