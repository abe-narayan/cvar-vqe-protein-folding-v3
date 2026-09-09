"""Regression tests for `s12/instrument.py` — the evaluation instrument.

WHY THIS FILE EXISTS
====================
Every RMSD in this project — every finding in `FINDINGS.md`, every sprint LEDGER row,
every incumbent number — is produced by a handful of functions in `s12/instrument.py`:
`kabsch_rmsd_batch`, `superpose_batch`, `pairwise_rmsd`, `medoid`, `coordinate_average`,
`pair_index`, `pair_dists`, `shipped_score`, `paired` and `write`.  Until now the module
had **no test of its own**.  `tests/test_amber_frame_invariance.py` imports it, but only
to obtain targets; it asserts nothing about the instrument's arithmetic.

That is the highest-leverage untested surface in the repository: a silent change to
`kabsch_rmsd_batch` would move every number in the project at once and nothing would fail.

THE DUPLICATE-UTILITY HAZARD THIS FILE PINS
===========================================
`s12/instrument.py` carries its OWN implementations of four functions that also exist in
`core.geometry`: `kabsch_rmsd_batch`, `ca_rmsd`, `pair_index`, `pair_dists`.  They are not
imports and not aliases — they are separate code that happens to agree today.  The
production pipeline uses `core.geometry`'s; the instrument that scores the pipeline uses
its own.  **If those two ever diverge, the project would be scoring one geometry with
another geometry's ruler and no existing test would notice.**  The cross-checks below
assert agreement to a stated tolerance, and are the reason the divergence cannot happen
quietly.  (They do NOT assert bit-equality: the two use different reduction orders — the
instrument reconstructs the RMSD from singular values, `core.geometry` from superposed
coordinates — so a few ulp of difference is expected and is measured rather than assumed.)

Everything here is pure arithmetic on synthetic arrays.  Nothing loads the window
universes, the distogram, ESM, OpenMM or any cache, so the file runs in about a second and
is safe to run while experiments hold the heavy locks.

Run:  python -m pytest tests/test_instrument.py -q
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from s12 import instrument as I                                       # noqa: E402

#: MEASURED, not assumed: on every case below the instrument and `core.geometry` agree
#: to EXACTLY 0.0, because both use the singular-value trace form of the residual rather
#: than forming the rotation.  The assertions are made at `== 0.0` where that holds, and
#: this tolerance is used only where a different reduction order is genuinely involved.
CROSS_TOL = 1e-9

#: THE INSTRUMENT'S NUMERICAL FLOOR, measured here rather than wished away.
#: `kabsch_rmsd_batch` reconstructs the residual as ``|P|^2 + |T|^2 - 2*sum(s)``.  When the
#: two structures are the same, that is a difference of two quantities of order 1e3-1e4 Å^2
#: whose true difference is zero, so the cancellation leaves ~1e-11 Å^2 of rounding, and the
#: square root turns it into ~1e-7 Å.  `core.geometry` uses the identical construction and
#: has the identical floor.  It is 4 orders of magnitude below the 1e-3 Å the project
#: reports to and cannot affect a finding — but an "RMSD of a structure with itself" is
#: 1.3e-7, never 0.0, and a test written to expect an exact zero is the test that is wrong.
SELF_RMSD_FLOOR = 1e-6


# ------------------------------------------------------------------ helpers
def _rng(seed=0):
    return np.random.default_rng(seed)


def _rot(rng):
    """A proper rotation (det = +1), from a QR of a Gaussian matrix."""
    Q, R = np.linalg.qr(rng.standard_normal((3, 3)))
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1.0
    return Q


def _chain(rng, n=12, b=None):
    """A plausible CA trace (or batch of them): a random walk at ~3.8 Å spacing."""
    shape = (n, 3) if b is None else (b, n, 3)
    step = rng.standard_normal(shape)
    step /= np.linalg.norm(step, axis=-1, keepdims=True)
    return np.cumsum(step * 3.8, axis=-2)


# ================================================================== kabsch_rmsd_batch
def test_kabsch_rmsd_of_a_structure_with_itself_is_zero_to_the_floor():
    """See `SELF_RMSD_FLOOR`.  The self-RMSD is ~1e-7 Å, not 0.0, and that is a property of
    the trace form both implementations use, not a defect.  Pinned so that a future change
    to a coordinate-difference form (which would give ~1e-15) is noticed as a change."""
    W = _chain(_rng(1), n=14, b=5)
    for k in range(len(W)):
        got = I.kabsch_rmsd_batch(W[k][None], W[k])[0]
        assert 0.0 <= got < SELF_RMSD_FLOOR


def test_kabsch_rmsd_is_invariant_under_rigid_motion_of_either_argument():
    """The whole point of superposition: translating or rotating an input must not move
    the answer.  Asserted on BOTH arguments, because a bug in the centring of one of them
    is invisible if only the other is tested."""
    rng = _rng(2)
    W, T = _chain(rng, 13, 6), _chain(rng, 13)
    base = I.kabsch_rmsd_batch(W, T)
    R, t = _rot(rng), rng.standard_normal(3) * 17.0
    assert I.kabsch_rmsd_batch(W @ R.T + t, T) == pytest.approx(base, abs=1e-9)
    assert I.kabsch_rmsd_batch(W, T @ R.T + t) == pytest.approx(base, abs=1e-9)


def test_kabsch_rmsd_is_symmetric_in_its_two_arguments():
    rng = _rng(3)
    A, B = _chain(rng, 11), _chain(rng, 11)
    assert I.kabsch_rmsd_batch(A[None], B)[0] == pytest.approx(
        I.kabsch_rmsd_batch(B[None], A)[0], abs=1e-9)


def test_kabsch_rmsd_refuses_the_reflection():
    """A mirror image is NOT the same structure.

    The determinant correction in the SVD is the single line that separates a proper
    rotation from an improper one.  Delete it and a mirrored chain superposes onto its
    original with RMSD 0, which would silently halve the reported error of every
    left-handed decoy in the project.  So the test asserts the mirror does NOT fit.
    """
    rng = _rng(4)
    T = _chain(rng, 12)
    M = T * np.array([1.0, 1.0, -1.0])          # reflection through the xy plane
    assert I.kabsch_rmsd_batch(M[None], T)[0] > 1.0


def test_kabsch_rmsd_matches_an_independent_closed_form():
    """Checked against a direct Kabsch written here from the textbook recipe (centre,
    covariance, SVD, det fix, rotate, then the plain RMSD of the moved coordinates).
    That reference shares no code with the instrument: the instrument reconstructs the
    RMSD from the singular values without ever forming the rotation."""
    rng = _rng(5)
    W, T = _chain(rng, 15, 7), _chain(rng, 15)

    def reference(P, Q):
        P = P - P.mean(0)
        Q = Q - Q.mean(0)
        U, S, Vt = np.linalg.svd(P.T @ Q)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        return float(np.sqrt((((R @ P.T).T - Q) ** 2).sum() / len(P)))

    got = I.kabsch_rmsd_batch(W, T)
    want = np.array([reference(w, T) for w in W])
    assert got == pytest.approx(want, abs=1e-9)


def test_kabsch_rmsd_scales_linearly_with_the_structure():
    rng = _rng(6)
    W, T = _chain(rng, 10, 4), _chain(rng, 10)
    assert I.kabsch_rmsd_batch(3.0 * W, 3.0 * T) == pytest.approx(
        3.0 * I.kabsch_rmsd_batch(W, T), rel=1e-9)


def test_kabsch_rmsd_batch_equals_looping_one_at_a_time():
    """Batching is an optimisation; it must not be a different answer."""
    rng = _rng(7)
    W, T = _chain(rng, 9, 11), _chain(rng, 9)
    batched = I.kabsch_rmsd_batch(W, T)
    looped = np.array([I.kabsch_rmsd_batch(w[None], T)[0] for w in W])
    assert batched == pytest.approx(looped, abs=1e-12)


def test_ca_rmsd_is_the_single_structure_face_of_the_batch_routine():
    rng = _rng(8)
    a, b = _chain(rng, 13), _chain(rng, 13)
    assert I.ca_rmsd(a, b) == I.kabsch_rmsd_batch(a[None], b)[0]


# ================================================================== cross-check vs core
def test_instrument_kabsch_agrees_with_core_geometry():
    """THE DIVERGENCE ALARM.  See the module docstring: the instrument scores what
    `core.geometry` builds, using its own copy of the same maths.  If this ever fails,
    every RMSD in the project is being measured with a different ruler than the one the
    pipeline optimises against, and no other test in the repository would say so."""
    from core import geometry as geo
    rng = _rng(9)
    W, T = _chain(rng, 14, 8), _chain(rng, 14)
    # Asserted at EXACT equality, not a tolerance: both reconstruct the residual from the
    # singular values in the same order, so the difference is measured at exactly 0.0 and
    # a tolerance here would only hide the first ulp of a real divergence.
    assert np.array_equal(I.kabsch_rmsd_batch(W, T),
                          np.asarray(geo.kabsch_rmsd_batch(W, T), float))
    # ...including on near-degenerate input, where the cancellation is worst.
    assert np.array_equal(I.kabsch_rmsd_batch(W, W[0]),
                          np.asarray(geo.kabsch_rmsd_batch(W, W[0]), float))


def test_instrument_pair_index_and_pair_dists_agree_with_core_geometry():
    from core import geometry as geo
    rng = _rng(10)
    for n in (8, 12, 17):
        for sep in (1, 2, 3):
            ii, jj = I.pair_index(n, sep)
            gi, gj = geo.pair_index(n, sep)
            assert np.array_equal(ii, gi) and np.array_equal(jj, gj)
            W = _chain(rng, n, 5)
            assert I.pair_dists(W, ii, jj) == pytest.approx(
                np.asarray(geo.pair_dists(W, gi, gj), float), abs=CROSS_TOL)


# ================================================================== superpose_batch
def test_superpose_batch_realises_the_rmsd_it_reports():
    """`superpose_batch` and `kabsch_rmsd_batch` are two functions that must describe the
    same superposition; the plain RMSD of the moved coordinates has to reproduce the
    value the other routine returns without moving anything."""
    rng = _rng(11)
    W, T = _chain(rng, 12, 6), _chain(rng, 12)
    moved = I.superpose_batch(W, T)
    direct = np.sqrt(((moved - T) ** 2).sum((1, 2)) / W.shape[1])
    assert direct == pytest.approx(I.kabsch_rmsd_batch(W, T), abs=1e-9)


def test_superpose_batch_is_a_rigid_motion():
    """Superposition may not deform.  Every internal distance survives it exactly."""
    rng = _rng(12)
    W, T = _chain(rng, 11, 4), _chain(rng, 11)
    ii, jj = I.pair_index(11, 1)
    before = I.pair_dists(W, ii, jj)
    after = I.pair_dists(I.superpose_batch(W, T), ii, jj)
    assert after == pytest.approx(before, abs=1e-9)


def test_superpose_batch_puts_the_centroid_on_the_targets():
    rng = _rng(13)
    W, T = _chain(rng, 10, 5), _chain(rng, 10)
    assert I.superpose_batch(W, T).mean(1) == pytest.approx(
        np.tile(T.mean(0), (5, 1)), abs=1e-9)


# ================================================================== pairwise / medoid
def test_pairwise_rmsd_is_symmetric_with_a_zero_diagonal():
    rng = _rng(14)
    W = _chain(rng, 12, 7)
    P = I.pairwise_rmsd(W)
    assert P.shape == (7, 7)
    assert np.all(np.diag(P) < SELF_RMSD_FLOOR)          # the floor, not exactly 0.0
    assert P == pytest.approx(P.T, abs=1e-9)


def test_pairwise_rmsd_is_not_symmetrised_the_way_core_geometry_symmetrises_it():
    """A real behavioural difference between the two implementations, pinned rather than
    tidied: `core.geometry.pairwise_ca_rmsd` returns ``0.5*(M + M.T)``, the instrument
    returns the raw matrix.  Both are symmetric to ~1e-9 here, so the choice does not move
    a number today — but a medoid is an argmin over ROW MEANS, and averaging the transpose
    in changes those means in the last bits, which is exactly the kind of difference that
    decides a tie.  Anyone unifying the two must know which one the findings were measured
    with: the instrument's, unsymmetrised."""
    from core import geometry as geo
    W = _chain(_rng(30), n=11, b=6)
    P_i, P_c = I.pairwise_rmsd(W), np.asarray(geo.pairwise_ca_rmsd(W), float)
    assert P_i == pytest.approx(P_c, abs=1e-9)
    assert P_c == pytest.approx(0.5 * (P_c + P_c.T), abs=0.0)


def test_pairwise_rmsd_entries_are_the_batch_routines_entries():
    rng = _rng(15)
    W = _chain(rng, 9, 6)
    P = I.pairwise_rmsd(W)
    for a in range(6):
        assert P[a] == pytest.approx(I.kabsch_rmsd_batch(W, W[a]), abs=1e-12)


def test_medoid_is_the_argmin_of_the_mean_row():
    rng = _rng(16)
    P = np.abs(rng.standard_normal((9, 9)))
    P = (P + P.T) / 2.0
    np.fill_diagonal(P, 0.0)
    assert I.medoid(P) == int(np.argmin(P.mean(1)))


def test_medoid_tie_break_is_the_lowest_index_and_that_is_a_documented_hazard():
    """`medoid` is `np.argmin`, so a tie resolves to the LOWEST INDEX — which is the
    order the pool arrived in.

    The project has already been bitten by this once: an `np.argmin` over a tied
    native-free signal silently read the ORACLE sort order of the candidate list and
    invented a 1.386 Å "winner" that did not exist.  This test does not call the behaviour
    wrong — deterministic tie-breaking is the right default for reproducibility — it PINS
    it, so that any caller ranking on a signal with large tie sets knows it is inheriting
    the candidate order and must average the outcome over the tied argmin set instead.
    """
    P = np.ones((5, 5))
    np.fill_diagonal(P, 0.0)                    # every row mean identical
    assert I.medoid(P) == 0


def test_coordinate_average_agrees_with_its_own_documented_recipe():
    """S8-11's operator, spelled out: pick the medoid, superpose everything on it, mean."""
    rng = _rng(17)
    W = _chain(rng, 13, 8)
    C, b = I.coordinate_average(W)
    assert b == I.medoid(I.pairwise_rmsd(W))
    assert C == pytest.approx(I.superpose_batch(W, W[b]).mean(0), abs=1e-12)
    assert C.shape == (13, 3)


def test_coordinate_average_of_one_structure_is_that_structure():
    rng = _rng(18)
    W = _chain(rng, 10, 1)
    C, b = I.coordinate_average(W)
    assert b == 0
    assert C == pytest.approx(W[0], abs=1e-9)


def test_coordinate_average_of_rigid_copies_is_that_structure():
    """Averaging identical geometry in different frames must return the geometry, not a
    collapsed blob.  This is the control for the known contraction effect: averaging
    genuinely different candidates shortens the backbone by ~26%, so a test that only
    checked "the average is shorter" would pass on a broken implementation too."""
    rng = _rng(19)
    T = _chain(rng, 12)
    W = np.stack([T @ _rot(rng).T + rng.standard_normal(3) * 9.0 for _ in range(6)])
    C, _ = I.coordinate_average(W)
    assert I.ca_rmsd(C, T) == pytest.approx(0.0, abs=1e-8)


def test_coordinate_average_accepts_a_precomputed_distance_matrix():
    rng = _rng(20)
    W = _chain(rng, 11, 5)
    P = I.pairwise_rmsd(W)
    assert I.coordinate_average(W, P)[0] == pytest.approx(
        I.coordinate_average(W)[0], abs=1e-12)


# ================================================================== pair geometry
def test_pair_index_respects_min_sep_and_covers_every_admissible_pair():
    for n in (5, 9, 16):
        for sep in (1, 2, 4):
            i, j = I.pair_index(n, sep)
            assert np.all(j - i >= sep)
            assert len(i) == sum(1 for a in range(n) for b in range(n) if b - a >= sep)
            assert len(set(zip(i.tolist(), j.tolist()))) == len(i)


def test_pair_dists_matches_a_brute_force_loop_and_broadcasts():
    rng = _rng(21)
    n = 10
    i, j = I.pair_index(n, 2)
    W = _chain(rng, n, 4)
    D = I.pair_dists(W, i, j)
    assert D.shape == (4, len(i))
    for b in range(4):
        for k in range(len(i)):
            assert D[b, k] == pytest.approx(
                float(np.linalg.norm(W[b, i[k]] - W[b, j[k]])), abs=1e-12)
    assert I.pair_dists(W[0], i, j) == pytest.approx(D[0], abs=1e-12)   # unbatched


def test_pair_dists_is_invariant_under_rigid_motion():
    rng = _rng(22)
    i, j = I.pair_index(12, 2)
    W = _chain(rng, 12, 3)
    R, t = _rot(rng), rng.standard_normal(3) * 5.0
    assert I.pair_dists(W @ R.T + t, i, j) == pytest.approx(
        I.pair_dists(W, i, j), abs=1e-9)


# ================================================================== shipped_score
def _fake_dg(npairs=6, lo=2.0, step=0.05, ngrid=760):
    """A distogram-shaped risk table with a known argmin per pair, so the lookup can be
    checked without loading a model."""
    grid = lo + step * np.arange(ngrid)
    risk = np.zeros((npairs, ngrid))
    centres = np.linspace(5.0, 20.0, npairs)
    for p in range(npairs):
        risk[p] = np.abs(grid - centres[p])
    return {"grid": grid, "risk": risk}, centres


def test_shipped_score_is_the_mean_risk_over_pairs_at_the_looked_up_bin():
    dg, centres = _fake_dg()
    D = np.tile(centres, (3, 1))
    D[1] += 1.0
    D[2] -= 2.0
    got = I.shipped_score(dg, D)
    assert got[0] == pytest.approx(0.0, abs=0.05)      # at the minima
    assert got[1] == pytest.approx(1.0, abs=0.05)
    assert got[2] == pytest.approx(2.0, abs=0.05)
    assert got[0] < got[1] < got[2]                    # lower is better, and it orders


def test_shipped_score_clamps_out_of_grid_distances_instead_of_wrapping():
    """Distances outside 2–40 Å do occur — a stretched decoy, a mis-built chain — and the
    lookup index is `clip`ped.  A wrap (a negative index) would silently score a 60 Å pair
    with the risk of a 40 Å one at the far END of the table, which is a different number.
    Both ends are asserted, because only the low end can go negative."""
    dg, _ = _fake_dg()
    ngrid = len(dg["grid"])
    below = I.shipped_score(dg, np.full((1, 6), -50.0))
    at_lo = I.shipped_score(dg, np.full((1, 6), float(dg["grid"][0])))
    above = I.shipped_score(dg, np.full((1, 6), 1e6))
    at_hi = I.shipped_score(dg, np.full((1, 6), float(dg["grid"][ngrid - 1])))
    assert below == pytest.approx(at_lo, abs=1e-12)
    assert above == pytest.approx(at_hi, abs=1e-12)


def test_shipped_score_returns_one_number_per_candidate():
    dg, centres = _fake_dg()
    assert I.shipped_score(dg, np.tile(centres, (17, 1))).shape == (17,)


# ================================================================== statistics
def test_paired_reports_the_arithmetic_it_claims_to():
    rng = _rng(23)
    a = rng.standard_normal(60) + 3.0
    b = rng.standard_normal(60) + 3.2
    out = I.paired(a, b, n_boot=200, seed=1)
    d = a - b
    assert out["n"] == 60
    assert out["mean_a"] == pytest.approx(float(a.mean()), abs=1e-12)
    assert out["mean_b"] == pytest.approx(float(b.mean()), abs=1e-12)
    assert out["mean_diff"] == pytest.approx(float(d.mean()), abs=1e-12)
    assert out["median_diff"] == pytest.approx(float(np.median(d)), abs=1e-12)
    assert out["se"] == pytest.approx(float(d.std(ddof=1) / math.sqrt(60)), abs=1e-12)
    assert out["n_better"] == int((d < 0).sum())
    assert out["n_worse"] == int((d > 0).sum())


def test_paired_se_is_the_one_the_mde_rule_consumes():
    """The project's hard rule is `MDE = 2.8016 x SE`, PER COMPARISON.  That makes the
    `se` field load-bearing: it is the ddof=1 standard error OF THE PAIRED DIFFERENCE, not
    of either arm, and not the pooled two-sample SE.  Using an arm's SE here would misprice
    every MDE in every LEDGER by the amount of the pairing correlation."""
    rng = _rng(24)
    common = rng.standard_normal(40) * 5.0          # a large shared component
    a = common + rng.standard_normal(40) * 0.1
    b = common + rng.standard_normal(40) * 0.1
    out = I.paired(a, b, n_boot=100, seed=0)
    arm_se = float(a.std(ddof=1) / math.sqrt(40))
    assert out["se"] < arm_se / 5.0                 # pairing must actually pay
    assert out["se"] == pytest.approx(
        float((a - b).std(ddof=1) / math.sqrt(40)), abs=1e-12)


def test_paired_is_antisymmetric_and_negative_means_a_is_better():
    rng = _rng(25)
    a = rng.standard_normal(30)
    b = a + 0.5                                     # a is uniformly lower == better
    out = I.paired(a, b, n_boot=100, seed=0)
    assert out["mean_diff"] < 0
    assert out["n_better"] == 30 and out["n_worse"] == 0
    rev = I.paired(b, a, n_boot=100, seed=0)
    assert rev["mean_diff"] == pytest.approx(-out["mean_diff"], abs=1e-12)


def test_paired_is_deterministic_given_the_seed():
    rng = _rng(26)
    a, b = rng.standard_normal(50), rng.standard_normal(50)
    assert I.paired(a, b, n_boot=300, seed=7) == I.paired(a, b, n_boot=300, seed=7)
    assert I.paired(a, b, n_boot=300, seed=7)["ci95"] != \
        I.paired(a, b, n_boot=300, seed=8)["ci95"]


def test_paired_ci_brackets_the_mean_and_narrows_with_n():
    rng = _rng(27)
    small = I.paired(rng.standard_normal(25), rng.standard_normal(25),
                     n_boot=2000, seed=3)
    big = I.paired(rng.standard_normal(2000), rng.standard_normal(2000),
                   n_boot=2000, seed=3)
    for out in (small, big):
        lo, hi = out["ci95"]
        assert lo <= out["mean_diff"] <= hi
    assert (big["ci95"][1] - big["ci95"][0]) < (small["ci95"][1] - small["ci95"][0])


def test_paired_reports_the_concentration_diagnostics_the_ledgers_quote():
    """The median-vs-mean gap and the drop-top means are the project's free early warning
    for a result carried by a handful of targets.  They must be the mean over the SAME
    difference vector with the largest gains removed — not a re-bootstrap, not a trim of
    both tails."""
    rng = _rng(28)
    a = rng.standard_normal(50)
    b = a.copy()
    b[:5] += 10.0                                   # five targets carry everything
    out = I.paired(a, b, n_boot=100, seed=0, names=[f"T{k}" for k in range(50)])
    d = a - b
    order = np.argsort(d)
    assert out["drop_top10_mean_diff"] == pytest.approx(float(d[order[10:]].mean()),
                                                        abs=1e-12)
    assert out["drop_top20_mean_diff"] == pytest.approx(float(d[order[20:]].mean()),
                                                        abs=1e-12)
    assert out["top10_share"] == pytest.approx(float(d[order[:10]].sum() / d.sum()),
                                               abs=1e-12)
    assert out["drop_top10_mean_diff"] > out["mean_diff"]      # the effect was the top 5
    assert len(out["top10_targets"]) == 10
    assert out["top10_targets"][0][0] in {f"T{k}" for k in range(5)}


def test_paired_per_fold_means_partition_the_targets():
    rng = _rng(29)
    a, b = rng.standard_normal(50), rng.standard_normal(50)
    folds = rng.integers(0, 5, 50)
    out = I.paired(a, b, n_boot=50, seed=0, folds=folds)
    d = a - b
    assert set(out["per_fold"]) == set(int(f) for f in np.unique(folds))
    for f, m in out["per_fold"].items():
        assert m == pytest.approx(float(d[folds == f].mean()), abs=1e-12)


def test_summary_reports_the_fields_the_findings_quote():
    x = np.array([1.0, 1.5, 2.5, 3.0, 9.0])
    s = I.summary(x)
    assert s["n"] == 5
    assert s["mean"] == pytest.approx(3.4)
    assert s["median"] == pytest.approx(2.5)
    assert s["min"] == 1.0 and s["max"] == 9.0
    assert s["sd"] == pytest.approx(float(x.std()))          # population sd, ddof=0
    assert s["frac_under_2.0"] == pytest.approx(0.4)


# ================================================================== the writer hazard
def test_write_flags_a_partial_result_as_incomplete(tmp_path, monkeypatch):
    """The documented hazard in `instrument.write`: it has NO CONFIG KEY, so a partial run
    overwrites a complete one of the same name.  The `complete` flag is the only thing that
    lets a reader tell them apart, and it must be driven by the ROW COUNT against
    `n_expected` — the project's provenance rule is that `complete` requires the full key
    set and never a filename."""
    monkeypatch.setattr(I, "RESULTS", str(tmp_path))
    part = json.load(open(I.write("x", {"per_target": [1, 2, 3]}, n_expected=126)))
    assert part["complete"] is False and part["n_rows"] == 3 and part["n_expected"] == 126
    full = json.load(open(I.write("y", {"per_target": list(range(126))},
                                  n_expected=126)))
    assert full["complete"] is True and full["n_rows"] == 126


def test_write_reads_rows_under_either_documented_key(tmp_path, monkeypatch):
    monkeypatch.setattr(I, "RESULTS", str(tmp_path))
    assert json.load(open(I.write("r", {"rows": [1, 2]}, n_expected=2)))["complete"]


def test_write_without_n_expected_makes_no_completeness_claim(tmp_path, monkeypatch):
    """Silence is better than a false `complete: true`."""
    monkeypatch.setattr(I, "RESULTS", str(tmp_path))
    assert "complete" not in json.load(open(I.write("z", {"per_target": [1]})))


def test_write_serialises_numpy_without_raising(tmp_path, monkeypatch):
    monkeypatch.setattr(I, "RESULTS", str(tmp_path))
    obj = json.load(open(I.write("n", {"per_target": [], "a": np.arange(3),
                                       "b": np.float64(1.5)}, n_expected=0)))
    assert obj["a"] == [0, 1, 2] and obj["b"] == 1.5


def test_write_appends_the_json_suffix_exactly_once(tmp_path, monkeypatch):
    monkeypatch.setattr(I, "RESULTS", str(tmp_path))
    assert I.write("a", {}).endswith("a.json")
    assert I.write("b.json", {}).endswith("b.json")
    assert not I.write("b.json", {}).endswith("b.json.json")


# ================================================================== constants
def test_the_pinned_constants_are_the_ones_every_finding_assumes():
    """K, M and BAND appear by value in the findings text.  If one of them moves, the
    numbers in `FINDINGS.md` stop referring to the experiment that produced them."""
    assert I.K == 500                     # the retrieval pool
    assert I.M == 75                      # the top-M that survive the filter
    assert I.BAND == 1.5                  # the near-native band the recall is defined on
    assert len(I.FAIL18) == 18 == len(set(I.FAIL18))


def test_the_alphabet_is_the_canonical_one_and_decode_inverts_core_encode():
    from core import data as D
    assert I.ALPHABET == D.ALPHABET
    seq = "ACDEFGHIKLMNPQRSTVWY"
    assert I.decode(np.asarray(D.encode(seq))) == seq
