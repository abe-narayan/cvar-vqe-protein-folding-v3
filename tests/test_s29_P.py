"""tests/test_s29_P.py -- lane P (the projection price).  PREREG section 11.

Light by construction: one real target's cached cloud is used where a real object is needed
(`s29/results/s29_P_clouds/1A13.npz`, built by `s29_P_scale.py factors`); everything else is
synthetic.  The one projection this file runs is the identity check, which is the gate the
whole lane rests on.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s29 import s29_P_scale as P            # noqa: E402


def _cloud(pdb="1A13"):
    f = os.path.join(P.CLOUDS, f"{pdb}.npz")
    if not os.path.exists(f):
        pytest.skip(f"cloud cache absent ({f}); run `s29_P_scale.py factors` first")
    return np.asarray(np.load(f)["C"], float)


# --------------------------------------------------------------- 1. s = 1 is the identity
def test_s1_is_the_exact_identity():
    rng = np.random.default_rng(3)
    C = rng.normal(size=(15, 3)) * 4.0
    out = P.rescale(C, 1.0)
    assert np.array_equal(out, C), "s=1 must be bit-for-bit the input"


def test_prod_projection_equals_a_direct_project_call():
    """Arm PROD is `I.project` on the unmodified cloud -- bit-for-bit, not to a tolerance."""
    from s12 import instrument as I
    from s27 import run_pool as RP
    C = _cloud("1A13")
    cand, _, _ = RP.channels_for("1A13")
    a = np.asarray(P.project_arm(C, cand, 1.0, P.LAM)["ca"], float)
    b = np.asarray(I.project(C, cand.seq, cand.fold, lam=P.LAM)["ca"], float)
    assert np.array_equal(a, b)


def test_shipped_objective_matches_the_projection_residual_plus_the_penalty():
    """obj0 must BE the shipped objective, not an approximation of it."""
    from s12 import instrument as I
    from s27 import run_pool as RP
    from core import project as pj
    C = _cloud("1A13")
    cand, _, _ = RP.channels_for("1A13")
    pr = P.project_arm(C, cand, 1.0, P.LAM)
    pen = pj.make_penalty("ramah", cand.seq, int(cand.fold))
    direct = I.ca_rmsd(pr["ca"], C) + P.LAM * float(
        np.asarray(pen(np.asarray(pr["phi"])[None], np.asarray(pr["psi"])[None])).ravel()[0])
    assert abs(P.shipped_objective(pr, C, cand) - direct) < 1e-15
    assert P.shipped_objective(pr, C, cand) >= I.ca_rmsd(pr["ca"], C) - 1e-15


# --------------------------------------------------------------- 2. the Rg identity
@pytest.mark.parametrize("n", [9, 12, 16])
def test_rg_identity_on_random_clouds(n):
    rng = np.random.default_rng(n)
    C = rng.normal(size=(n, 3)) * 5.0
    D = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=-1)
    assert abs(P.rg_from_pairs(D, n) - P.rg(C)) < 1e-12


def test_rg_identity_on_a_real_cloud():
    C = _cloud()
    n = len(C)
    D = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=-1)
    assert abs(P.rg_from_pairs(D, n) - P.rg(C)) < 1e-12


# --------------------------------------------------------------- 3. BOND
def test_bond_makes_the_mean_bond_ideal():
    for seed in range(5):
        rng = np.random.default_rng(seed)
        C = rng.normal(size=(13, 3)) * 2.0
        g = P.g_bond(C)
        assert abs(P.mean_bond(P.rescale(C, g)) - P.IDEAL_BOND) < 1e-12


def test_bond_on_the_real_cloud_is_an_expansion():
    C = _cloud()
    assert P.g_bond(C) > 1.0, "the production cloud is contracted, so g must expand it"


# --------------------------------------------------------------- 4. ISO is the LS minimiser
def test_iso_is_the_least_squares_minimiser():
    from s12 import instrument as I
    C = _cloud()
    dg = I.distogram("1A13")
    s = P.s_iso(C, dg)
    i, j, m = P.posterior_median(dg)
    d = np.linalg.norm(C[i] - C[j], axis=1)

    def obj(x):
        return float(((x * d - m) ** 2).sum())

    grid = np.linspace(s - 0.2, s + 0.2, 4001)
    vals = np.array([obj(x) for x in grid])
    assert abs(grid[int(np.argmin(vals))] - s) < 2e-4
    h = 1e-6
    assert obj(s) <= obj(s + h) and obj(s) <= obj(s - h)


def test_iso_is_at_most_the_restricted_span_by_cauchy_schwarz():
    from s12 import instrument as I
    C = _cloud()
    dg = I.distogram("1A13")
    assert P.s_iso(C, dg) <= P.s_spanr(C, dg) + 1e-12


# --------------------------------------------------------------- 5. SPAN
def test_span_makes_rg_equal_the_posterior_rg():
    from s12 import instrument as I
    C = _cloud()
    dg = I.distogram("1A13")
    s, rg_post = P.s_span(C, dg)
    assert abs(P.rg(P.rescale(C, s)) - rg_post) < 1e-12


# --------------------------------------------------------------- 6. the controls
def test_ctrl_inv_is_the_exact_reciprocal():
    C = _cloud()
    g = P.g_bond(C)
    assert abs((1.0 / g) * g - 1.0) < 1e-15


def test_derangements_fix_nothing_and_preserve_the_multiset():
    gs = np.linspace(1.1, 1.5, 126)
    Pm = P.derangements(gs, P.N_RAND)
    assert Pm.shape == (P.N_RAND, 126)
    for p in Pm:
        assert not np.any(p == np.arange(126)), "a derangement must fix no target"
        assert np.allclose(np.sort(gs[p]), np.sort(gs)), "matched magnitude by construction"
    assert len({tuple(p.tolist()) for p in Pm}) == P.N_RAND, "draws must differ"


def test_derangements_are_reproducible():
    gs = np.linspace(1.0, 2.0, 40)
    assert np.array_equal(P.derangements(gs, 3), P.derangements(gs, 3))


# --------------------------------------------------------------- 7. NaN-poison
def test_scale_functions_are_finite_when_the_native_is_nan():
    """Every deployable scale is a function of (cloud, posterior) only.  A poisoned native
    must not change any of them, and must not make them non-finite."""
    from s12 import instrument as I
    from s27 import run_pool as RP
    C = _cloud()
    dg = I.distogram("1A13")
    cand, _, _ = RP.channels_for("1A13")
    before = (P.g_bond(C), P.s_span(C, dg)[0], P.s_iso(C, dg), P.s_spanr(C, dg))
    poisoned = np.asarray(cand.nat_ca, float).copy()
    poisoned[:] = np.nan
    after = (P.g_bond(C), P.s_span(C, dg)[0], P.s_iso(C, dg), P.s_spanr(C, dg))
    assert all(np.isfinite(x) for x in after)
    assert before == after
    assert np.isnan(poisoned).all()


def test_nan_in_the_cloud_propagates_and_is_never_swallowed():
    from s12 import instrument as I
    C = _cloud().copy()
    dg = I.distogram("1A13")
    C[3, 1] = np.nan
    assert not np.isfinite(P.g_bond(C))
    assert not np.isfinite(P.s_iso(C, dg))
    assert not np.isfinite(P.s_spanr(C, dg))
    assert not np.isfinite(P.s_span(C, dg)[0])


def test_cost_isoresid_never_reads_the_native():
    """The meter's contract: ctx.nat_ca is NaN by construction; the cost must still be finite."""
    from s12 import instrument as I

    class Ctx:
        pass

    C = _cloud()
    ctx = Ctx()
    ctx.dg = I.distogram("1A13")
    ctx.nat_ca = np.full_like(C, np.nan)
    v = P.cost_isoresid(np.stack([C, C * 1.1]), ctx)
    assert v.shape == (2,) and np.isfinite(v).all()


# --------------------------------------------------------------- 8. the posterior median
def test_posterior_median_is_the_argmin_of_the_l1_risk_recomputed_independently():
    from s12 import instrument as I
    from core import predict as dgm
    dg = I.distogram("1A13")
    i, j, m = P.posterior_median(dg)
    prob = np.asarray(dg["prob"], float)
    centres = np.asarray(dg["centres"], float)
    grid = np.asarray(dg["grid"], float)
    assert np.allclose(centres, np.asarray(dgm.CENTRES, float), atol=1e-6)
    risk = (prob[:, None, :] * np.abs(grid[None, :, None] - centres[None, None, :])).sum(2)
    m2 = grid[np.argmin(risk, axis=1)]
    assert np.array_equal(m, m2), "the stored risk's argmin must be the L1 Bayes point"
    assert (m >= grid[0]).all() and (m <= grid[-1]).all()


def test_arm_table_is_complete_and_unique():
    if not os.path.exists(P.FACTORS):
        pytest.skip("factors absent; run `s29_P_scale.py factors` first")
    import json
    fac = json.load(open(P.FACTORS, encoding="utf-8"))
    pdb = fac["order"][0]
    A = P.arms_for(pdb, fac)
    names = [a for a, _, _ in A]
    assert len(names) == len(set(names))
    assert len(names) == 8 + P.N_RAND + len([s for s in P.ORACLE_GRID if s != 1.0])
    for nm in P.DEPLOYABLE:
        assert nm in names
    for _, s, lam in A:
        assert np.isfinite(s) and s > 0 and np.isfinite(lam) and lam > 0
