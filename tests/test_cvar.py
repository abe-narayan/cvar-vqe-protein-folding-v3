"""Definitional regression tests for the CVaR machinery in `core.quantum`.

WHAT THIS FILE ADDS THAT `tests/test_quantum.py` DOES NOT
=========================================================
`tests/test_quantum.py` is an EQUIVALENCE suite: it proves that the consolidated
`core.quantum` reproduces `qansatz` / `vqe` / `foldvqe` bit-for-bit, and it pins the one
gradient-baseline defect that was found and fixed.  That is the right test for a
consolidation, and it is thorough — but it inherits its notion of "correct" from the
module it is comparing against.  **If the shipped CVaR and its reference were both wrong
in the same way, every one of those tests would still pass.**

So this file asserts the DEFINITION instead.  CVaR_alpha of a distribution is the mean of
its lower alpha tail, and that object has properties that hold no matter who implemented
it: it equals the mean at alpha=1, it never exceeds the mean, it is monotone in alpha, it
descends to the minimum as alpha shrinks, it is translation-equivariant and positively
homogeneous, and the four entry points in this module — `cvar`, `cvar_from_samples`,
`cvar_from_probs`, `cvar_from_distribution` — must all agree wherever their domains meet.
None of that is checked anywhere else in the suite.

The four entry points, and why there are four:

    cvar(e, a)                  -> (value, quantile, boolean tail mask) over SAMPLES
    cvar_from_samples(e, a)     -> value only; the objective the SPSA driver evaluates
    cvar_from_probs(e, p, a)    -> EXACT, over an explicit distribution, tail mass split
                                   fractionally at the boundary state
    cvar_from_distribution(...) -> the same quantity by an explicit Python loop

The sampled pair and the exact pair are NOT the same function: the sampled pair takes a
whole number of samples (``k = ceil(alpha*N)``) and the exact pair splits the boundary
state's probability mass.  They coincide exactly when ``alpha*N`` is an integer and the
weights are uniform, and the tests below assert agreement there and only there — asserting
it everywhere would be asserting a falsehood.

Also covered here for the first time: `alpha_schedule` (the anneal that exists to prevent
premature concentration — a failure this project has met repeatedly), `tail_indices`
(the selection every gradient weight is gathered through), `BestSeenTracker` (which
carries the answer a run actually returns), `all_bitstrings` and `n_parameters`.

Pure arithmetic on synthetic arrays: no circuit is simulated, no device is opened, nothing
is loaded from disk.  Runs in under a second.

Run:  python -m pytest tests/test_cvar.py -q
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

import core.quantum as Q                                              # noqa: E402


def _rng(seed=0):
    return np.random.default_rng(seed)


def _reference_cvar(e, alpha):
    """The definition, written out: sort ascending, average the lowest ceil(alpha*N).

    Deliberately naive and deliberately not the shipped code path — no partition, no mask,
    no tie rule — so that it is an independent statement of what the answer should be.
    """
    e = np.sort(np.asarray(e, float))
    k = max(1, int(math.ceil(alpha * e.size)))
    return float(e[:k].mean())


# ===================================================== the definition, on samples
@pytest.mark.parametrize("alpha", [0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1.0])
@pytest.mark.parametrize("n", [1, 2, 7, 100, 449, 1000])
def test_cvar_is_the_mean_of_the_lower_tail_by_definition(alpha, n):
    """Both sampled entry points, against a naive independent reference, across the
    `CVAR_SORT_CUTOFF` boundary (448) so that BOTH the sort branch and the partition
    branch are exercised at every alpha."""
    e = _rng(n).standard_normal(n) * 3.0 + 1.0
    want = _reference_cvar(e, alpha)
    assert Q.cvar(e, alpha)[0] == pytest.approx(want, rel=1e-12, abs=1e-12)
    assert Q.cvar_from_samples(e, alpha) == pytest.approx(want, rel=1e-12, abs=1e-12)


def test_cvar_at_alpha_one_is_exactly_the_mean():
    """The boundary that matters most: if the whole optimiser silently ran at alpha=1 it
    would be minimising the MEAN, which is a different algorithm wearing CVaR's name."""
    e = _rng(1).standard_normal(500)
    assert Q.cvar(e, 1.0)[0] == pytest.approx(float(e.mean()), rel=1e-12)
    assert Q.cvar_from_samples(e, 1.0) == pytest.approx(float(e.mean()), rel=1e-12)
    assert Q.cvar(e, 1.0)[2].all()                       # the tail is everything


def test_cvar_never_exceeds_the_mean_and_never_undercuts_the_minimum():
    rng = _rng(2)
    for _ in range(20):
        e = rng.standard_normal(rng.integers(3, 400))
        for a in (0.05, 0.2, 0.6, 1.0):
            v = Q.cvar(e, a)[0]
            assert float(e.min()) - 1e-12 <= v <= float(e.mean()) + 1e-12


def test_cvar_is_monotone_non_decreasing_in_alpha():
    """A wider tail can only pull the average up.  This is the property that makes
    `alpha_schedule`'s anneal meaningful; if it failed, tightening alpha would not
    concentrate on the good states."""
    e = _rng(3).standard_normal(600)
    alphas = [0.02, 0.05, 0.1, 0.25, 0.5, 0.8, 1.0]
    vals = [Q.cvar(e, a)[0] for a in alphas]
    assert all(b >= a - 1e-12 for a, b in zip(vals, vals[1:])), vals


def test_cvar_descends_to_the_minimum_as_the_tail_narrows():
    e = _rng(4).standard_normal(1000)
    assert Q.cvar(e, 1.0 / len(e))[0] == pytest.approx(float(e.min()), rel=1e-12)


def test_cvar_is_translation_equivariant_and_positively_homogeneous():
    """CVaR(E + c) = CVaR(E) + c and CVaR(sE) = s CVaR(E) for s > 0.

    These are the two properties that let the project rank-standardise between
    Hamiltonians at all: an energy offset or a unit change must not reorder anything.
    A `max` where a `min` belongs, or an absolute value in the tail selection, breaks
    exactly these and nothing else in the suite would show it.
    """
    e = _rng(5).standard_normal(300) * 7.0
    for a in (0.1, 0.3, 1.0):
        base = Q.cvar(e, a)[0]
        assert Q.cvar(e + 12.5, a)[0] == pytest.approx(base + 12.5, rel=1e-12)
        assert Q.cvar(3.0 * e, a)[0] == pytest.approx(3.0 * base, rel=1e-12)
        # negation flips the tail to the OTHER end -- CVaR is not an even function
        assert Q.cvar(-e, a)[0] == pytest.approx(-_upper_tail_mean(e, a), rel=1e-12)


def _upper_tail_mean(e, alpha):
    e = np.sort(np.asarray(e, float))[::-1]
    k = max(1, int(math.ceil(alpha * e.size)))
    return float(e[:k].mean())


def test_cvar_of_a_constant_is_that_constant():
    for a in (0.01, 0.5, 1.0):
        assert Q.cvar(np.full(97, -4.25), a)[0] == pytest.approx(-4.25, rel=1e-15)


def test_cvar_selects_the_tail_and_not_the_mean_on_a_hand_worked_case():
    """Worked by hand so the expected numbers are readable rather than derived:
    energies 0..9, alpha=0.3 -> k=3 -> tail {0,1,2} -> CVaR 1.0, quantile 2.0.
    The mean is 4.5, so a routine that returned the mean fails loudly."""
    e = np.arange(10.0)
    v, q, tail = Q.cvar(e, 0.3)
    assert v == pytest.approx(1.0) and q == pytest.approx(2.0)
    assert tail.tolist() == [True] * 3 + [False] * 7
    assert v != pytest.approx(float(e.mean()))


# ===================================================== tail_indices
@pytest.mark.parametrize("n", [10, 449, 900])
@pytest.mark.parametrize("alpha", [0.05, 0.2, 0.5, 1.0])
def test_tail_indices_returns_exactly_k_elements_and_they_are_the_k_lowest(n, alpha):
    e = _rng(n + 7).standard_normal(n)
    k, q, tail = Q.tail_indices(e, alpha)
    assert k == max(1, int(math.ceil(alpha * n)))
    assert int(tail.sum()) == k
    # every selected element is <= every rejected one: this IS "the lowest k"
    if k < n:
        assert e[tail].max() <= e[~tail].min() + 1e-15
        assert q == pytest.approx(float(e[tail].max()))


def test_tail_indices_breaks_ties_by_lowest_index_on_both_branches():
    """Heavy ties are the case where a selection can silently disagree with a stable sort.
    The documented contract is `order[:k]` of a STABLE argsort, i.e. lowest index first.
    Checked on both sides of `CVAR_SORT_CUTOFF`, because the two branches reach that
    answer by different routes (argsort vs argpartition plus an explicit tie fixup)."""
    for n in (100, 2000):                       # below and above the 448 cutoff
        e = np.zeros(n)
        e[:n // 2] = -1.0                       # half tied low, half tied high
        k, q, tail = Q.tail_indices(e, 0.75)
        want = np.zeros(n, bool)
        want[np.argsort(e, kind="stable")[:k]] = True
        assert np.array_equal(tail, want), n


def test_tail_indices_mask_is_what_cvar_averages():
    e = _rng(9).standard_normal(500)
    for a in (0.1, 0.4):
        v, q, tail = Q.cvar(e, a)
        k2, q2, tail2 = Q.tail_indices(e, a)
        assert np.array_equal(tail, tail2) and q == q2
        assert v == pytest.approx(float(e[tail].mean()), rel=1e-15)


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.0001, float("nan")])
def test_cvar_entry_points_reject_an_alpha_outside_the_half_open_unit_interval(bad):
    e = np.arange(10.0)
    for fn in (Q.cvar, Q.cvar_from_samples, Q.tail_indices):
        with pytest.raises(ValueError):
            fn(e, bad)
    with pytest.raises(ValueError):
        Q.cvar_from_probs(e, np.full(10, 0.1), bad)


def test_cvar_entry_points_reject_an_empty_sample():
    for fn in (Q.cvar, Q.cvar_from_samples, Q.tail_indices):
        with pytest.raises(ValueError):
            fn(np.array([]), 0.1)


def test_cvar_from_probs_rejects_zero_total_mass():
    with pytest.raises(ValueError):
        Q.cvar_from_probs(np.arange(4.0), np.zeros(4), 0.5)


# ===================================================== the exact, distributional form
def test_cvar_from_probs_and_from_distribution_are_the_same_quantity():
    """Two implementations of the exact CVaR — one vectorised, one an explicit loop —
    and they exist so that the vectorised one can be differentiated.  If they disagree,
    the finite-difference gradient tests are differentiating something other than the
    objective the run minimises."""
    rng = _rng(11)
    for _ in range(25):
        n = int(rng.integers(2, 40))
        e = rng.standard_normal(n) * 4.0
        p = rng.random(n) + 1e-3
        p /= p.sum()
        for a in (0.05, 0.17, 0.5, 0.99, 1.0):
            assert Q.cvar_from_probs(e, p, a)[0] == pytest.approx(
                Q.cvar_from_distribution(e, p, a), rel=1e-9, abs=1e-12)


def test_exact_cvar_matches_the_sampled_one_when_the_boundary_does_not_split():
    """The one place the sampled and exact forms MUST coincide: uniform weights and
    `alpha*N` a whole number, so no state's mass is cut.  Asserting agreement anywhere
    else would be asserting something false — see the module docstring."""
    rng = _rng(12)
    n = 20
    e = rng.standard_normal(n)
    p = np.full(n, 1.0 / n)
    for k in (1, 2, 5, 10, 20):
        a = k / n
        assert Q.cvar_from_probs(e, p, a)[0] == pytest.approx(
            Q.cvar(e, a)[0], rel=1e-9, abs=1e-12)


def test_exact_cvar_tail_mass_totals_alpha_and_never_exceeds_a_states_own_mass():
    """The returned mass vector is what makes the value differentiable in `probs`.  It has
    to be a genuine sub-distribution: non-negative, capped by each state's probability,
    and summing to exactly alpha."""
    rng = _rng(13)
    for _ in range(20):
        n = int(rng.integers(3, 30))
        e = rng.standard_normal(n)
        p = rng.random(n) + 1e-3
        p /= p.sum()
        for a in (0.1, 0.44, 1.0):
            v, q, mass = Q.cvar_from_probs(e, p, a)
            assert mass.min() >= -1e-15
            assert np.all(mass <= p + 1e-12)
            assert mass.sum() == pytest.approx(a, abs=1e-12)
            assert v == pytest.approx(float((mass * e).sum() / a), rel=1e-12)


def test_exact_cvar_puts_all_its_mass_on_the_lowest_states():
    """Nothing may be taken from a high-energy state while a lower one is not yet full."""
    e = np.array([0.0, 1.0, 2.0, 3.0])
    p = np.array([0.25, 0.25, 0.25, 0.25])
    v, q, mass = Q.cvar_from_probs(e, p, 0.5)
    assert mass == pytest.approx([0.25, 0.25, 0.0, 0.0], abs=1e-12)
    assert v == pytest.approx(0.5)                       # (0 + 1) / 2


def test_exact_cvar_splits_the_boundary_state_fractionally():
    """Half of the third state's mass, and no more — the property that separates the exact
    form from the sampled one."""
    e = np.array([0.0, 1.0, 2.0, 3.0])
    p = np.array([0.25, 0.25, 0.25, 0.25])
    v, q, mass = Q.cvar_from_probs(e, p, 0.625)
    assert mass == pytest.approx([0.25, 0.25, 0.125, 0.0], abs=1e-12)
    assert v == pytest.approx((0.0 * .25 + 1.0 * .25 + 2.0 * .125) / 0.625)


def test_exact_cvar_at_alpha_one_is_the_full_expectation():
    rng = _rng(14)
    e = rng.standard_normal(30)
    p = rng.random(30)
    p /= p.sum()
    assert Q.cvar_from_probs(e, p, 1.0)[0] == pytest.approx(float((p * e).sum()),
                                                            rel=1e-12)


def test_exact_cvar_normalises_unnormalised_weights():
    e = np.array([0.0, 1.0, 2.0, 3.0])
    p = np.array([1.0, 1.0, 1.0, 1.0])
    assert Q.cvar_from_probs(e, p, 0.5)[0] == pytest.approx(
        Q.cvar_from_probs(e, p / 4.0, 0.5)[0], rel=1e-12)


def test_exact_cvar_is_translation_equivariant():
    rng = _rng(15)
    e, p = rng.standard_normal(25), rng.random(25)
    p /= p.sum()
    for a in (0.2, 0.7):
        assert Q.cvar_from_probs(e + 9.0, p, a)[0] == pytest.approx(
            Q.cvar_from_probs(e, p, a)[0] + 9.0, rel=1e-11)


def test_exact_cvar_ignores_states_of_zero_probability():
    """An unreachable state must not drag the objective down.  This is the failure mode
    where a selector 'finds' a brilliant configuration the circuit never actually prepares."""
    e = np.array([-99.0, 0.0, 1.0, 2.0])
    p = np.array([0.0, 1 / 3, 1 / 3, 1 / 3])
    v, q, mass = Q.cvar_from_probs(e, p, 0.5)
    assert mass[0] == pytest.approx(0.0, abs=1e-15)
    assert v > 0.0


# ===================================================== alpha_schedule
def test_alpha_schedule_hits_both_endpoints_exactly():
    assert Q.alpha_schedule(0.0) == pytest.approx(0.5)
    assert Q.alpha_schedule(1.0) == pytest.approx(0.05)
    assert Q.alpha_schedule(0.0, 0.4, 0.02) == pytest.approx(0.4)
    assert Q.alpha_schedule(1.0, 0.4, 0.02) == pytest.approx(0.02)


def test_alpha_schedule_is_monotone_decreasing_and_geometric():
    vals = [Q.alpha_schedule(p / 20.0) for p in range(21)]
    assert all(b < a for a, b in zip(vals, vals[1:]))
    # geometric: equal progress steps give a constant ratio
    ratios = [b / a for a, b in zip(vals, vals[1:])]
    assert np.allclose(ratios, ratios[0], rtol=1e-12)


def test_alpha_schedule_clamps_progress_outside_zero_one():
    assert Q.alpha_schedule(-5.0) == pytest.approx(Q.alpha_schedule(0.0))
    assert Q.alpha_schedule(17.0) == pytest.approx(Q.alpha_schedule(1.0))


def test_alpha_schedule_stays_a_legal_cvar_level():
    """Its whole output range has to be acceptable to the functions that consume it —
    a schedule that emitted 0.0 or >1 would raise deep inside a run."""
    for p in np.linspace(-1, 2, 61):
        a = Q.alpha_schedule(float(p))
        assert 0.0 < a <= 1.0
        Q.cvar(np.arange(20.0), a)                        # must not raise


# ===================================================== BestSeenTracker
def test_best_seen_tracker_keeps_the_minimum_and_counts_every_lookup():
    """The tracker carries what a run actually returns.  Two things are asserted because
    each has been a real bug class: it must keep the LOWEST (not the last, not the
    highest), and `n_lookups` must count every offer including rejected ones — that count
    is the evidence used to show the search does not enumerate `2**n`."""
    t = Q.BestSeenTracker()
    assert t.best_energy == float("inf") and t.best_bitstring is None
    for bits, e in [("0000", 3.0), ("1010", -1.5), ("1111", 7.0), ("0110", -0.5)]:
        t.offer(bits, e)
    assert t.best_energy == -1.5 and t.best_bitstring == "1010"
    assert t.n_lookups == 4


def test_best_seen_tracker_keeps_the_first_of_equal_energies():
    t = Q.BestSeenTracker()
    t.offer("aaa", 1.0)
    t.offer("bbb", 1.0)
    assert t.best_bitstring == "aaa"          # strict `<`, so a tie does not displace


# ===================================================== register helpers
def test_all_bitstrings_enumerates_the_register_msb_first():
    b = Q.all_bitstrings(3)
    assert b.shape == (8, 3)
    assert b[0].tolist() == [0, 0, 0]
    assert b[1].tolist() == [0, 0, 1]         # MSB first
    assert b[7].tolist() == [1, 1, 1]
    assert len({tuple(r) for r in b.tolist()}) == 8
    assert (b @ (1 << np.arange(2, -1, -1))).tolist() == list(range(8))


def test_all_bitstrings_refuses_to_become_load_bearing():
    """It is a verification instrument, and the ceiling is what keeps it one."""
    with pytest.raises(ValueError):
        Q.all_bitstrings(21)


def test_n_parameters_is_one_rotation_per_wire_per_layer():
    for nq in (1, 4, 12):
        for L in (1, 3, 7):
            assert Q.n_parameters(nq, L) == nq * L


# ===================================================== cross-entry-point consistency
def test_the_two_sampled_entry_points_never_disagree():
    """`cvar` gathers through a boolean mask; `cvar_from_samples` sorts the partitioned
    values and never forms an index.  They are separate code with separate tie rules and
    they must still return the same float, including under heavy ties where the tie rule
    could matter."""
    rng = _rng(16)
    for n in (5, 448, 449, 3000):
        for e in (rng.standard_normal(n),
                  np.round(rng.standard_normal(n) * 2),      # heavy ties
                  np.zeros(n)):
            for a in (0.03, 0.25, 0.5, 1.0):
                assert Q.cvar(e, a)[0] == pytest.approx(
                    Q.cvar_from_samples(e, a), rel=1e-12, abs=1e-12), (n, a)


def test_cvar_ignores_the_order_its_samples_arrive_in():
    """A permutation of the same measurements is the same distribution.  If this failed,
    a run's answer would depend on the order the shots came back in."""
    rng = _rng(17)
    e = rng.standard_normal(777)
    perm = rng.permutation(len(e))
    for a in (0.07, 0.4, 1.0):
        assert Q.cvar(e[perm], a)[0] == pytest.approx(Q.cvar(e, a)[0], rel=1e-12)
        assert Q.cvar_from_samples(e[perm], a) == pytest.approx(
            Q.cvar_from_samples(e, a), rel=1e-12)


def test_cvar_survives_the_scale_raw_amber_actually_arrives_at():
    """Unrelaxed AMBER windows are clash-dominated: 57.5% above 1e4 kcal and a worst case
    of 7.1e18.  The pipeline rank-standardises before scoring for exactly that reason, but
    the CVaR routine must not itself overflow, return NaN, or silently reorder when a
    handful of such values are present — the tail is the LOW end and must be unaffected."""
    e = np.concatenate([_rng(18).standard_normal(200), np.full(20, 7.1e18)])
    v, q, tail = Q.cvar(e, 0.1)
    assert np.isfinite(v)
    assert not tail[200:].any()                           # no clash entered the tail
    assert v == pytest.approx(_reference_cvar(e, 0.1), rel=1e-12)
