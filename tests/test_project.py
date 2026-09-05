"""Equivalence, validity and convention pins for `core.project` -- STAGE 3b.

The stage was rewritten for speed.  Speed is not the deliverable; EQUIVALENCE is.  Every
test here answers one of:

  (a) does the SHIPPED projection reproduce the reference on real pipeline inputs -- and
      the bar is EXACTLY 0.0, not "within noise"?
  (b) for the fast arm, which is a different pipeline and hashed as one: by how much does
      it move the answer, and on how many targets?
  (c) is the emitted structure still legal geometry and still L-handed?
  (d) is a parameter that is science rather than tuning -- the number of starts, maxiter,
      the penalty, lam -- still exactly what it was?

WHY (a) IS 0.0 AND NOT A TOLERANCE.  The first version of this file allowed the optimised
projection to land "within optimiser path noise" of the reference.  An independent audit
showed that was the wrong bar: a forward map agreeing with the reference to 1e-13 A still
emits a DIFFERENT STRUCTURE on 126/126 real targets, median 0.031 A and worst 1.63 A, with
the reference's own gradient formula held fixed.  The projection is degenerate -- two
ideal-geometry torsion branches at near-equal objective distance -- and L-BFGS-B is started
from a fully extended or fully helical chain, far from any minimum, where its early steps
are large enough to carry a 1e-13 A difference into the other branch.  So for this stage,
"numerically negligible" and "returns the same structure" are simply not the same claim,
and only the strict bar detects the difference.

This project has been burned by pre-projection "gains" of 0.05-0.12 A that sat on chains
with 2.65-3.13 A bonds and cost +0.30 to +0.95 A once legal geometry was restored, so a
CA-RMSD is never accepted here without the bond, angle, clash and chirality table beside
it.

The heavy per-target tables over all 126 tuning targets live in `verify/project_*.json`
(``equiv``, ``exactness``, ``stability``, ``degeneracy``, ``iters``); the tests read them
AND re-run a small live subset, so a stale file cannot carry the suite.

Run: ``python -m pytest tests/test_project.py -q``
"""
import json
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import core                                                            # noqa: E402
from core import project as pj                                         # noqa: E402

#: the arms `verify/project_equiv.json` carries, in report order
ARMS = ("ref", "ex", "fd", "an")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EQUIV = os.path.join(ROOT, "verify", "project_equiv.json")

#: The instrument the optimised path must reproduce, from the 126-target tuning run.
#: `EXPECT_FIT` is the multi-start lam=0 projection, which is what `s9/final.py` emits and
#: what benchmark60's pre-registration records as `tuning_mean`.  `EXPECT_SYNTHESIS` is the
#: shipped `ramah@0.3` arm.  Both are the REFERENCE's numbers, reproduced here, not targets.
EXPECT_FIT = 3.204
EXPECT_SYNTHESIS = 3.2148
EXPECT_RMSD_AVG = 3.0483
EXPECT_PROJECTION_COST = 0.166
#: Ideal CA-CA bond of the builder's chain, and the natives' mean pseudo-angle.
EXPECT_BOND = 3.804
NATIVE_PSEUDOANGLE = 103.99

#: The DEFAULT projection must be bit-identical to the reference.  Not "close": 0.0.
#: An independent audit showed why the weaker bar was wrong -- the fast scan builder agrees
#: with the reference to 1e-13 A and still moves the emitted structure on 126/126 targets,
#: because the projection is degenerate and L-BFGS-B starts far from any minimum.
TOL_EXACT = 0.0
#: The `analytic` arm is a DIFFERENT PIPELINE, hashed as one.  These bound what it is
#: allowed to do to the aggregate, not to any single target -- it demonstrably changes
#: single targets by up to 1.9 A and that is reported, not tolerated away.
TOL_MEAN_RMSD = 0.02        # A, mean |change in CA-RMSD to native|


# =====================================================================================
# helpers
# =====================================================================================
@pytest.fixture(scope="module")
def inputs():
    p = pj.INPUTS_JSON
    if not os.path.exists(p):
        pytest.skip("verify/project_inputs.json absent; run `python -m core.project harvest`")
    with open(p) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def measured():
    if not os.path.exists(EQUIV):
        pytest.skip("verify/project_equiv.json absent; run `python -m core.project equiv`")
    with open(EQUIV) as fh:
        return json.load(fh)


def _run(mod, row, grad=None):
    C = np.asarray(row["C"], float)
    pen = mod.make_penalty("ramah", row["seq"], int(row["fold"]))
    kw = {"grad": grad} if grad and mod is pj else {}
    p = mod.lam_path(C, pen, (0.0, 0.3), maxiter=300, multi=True, **kw)
    return p[0.0], p[0.3]


# =====================================================================================
# 1. THE FORWARD MAP
# =====================================================================================
def test_scan_builder_reproduces_the_reference_builder():
    """The scan is a different SUMMATION ORDER for the same geometry, so the difference is
    accumulated rounding and must scale like the chain, not like an error."""
    geo = core.backend("geometry")
    rng = np.random.default_rng(11)
    worst = {}
    for n in (2, 3, 9, 13, 21, 40, 80):
        phi = rng.uniform(-math.pi, math.pi, (6, n))
        psi = rng.uniform(-math.pi, math.pi, (6, n))
        ref = geo.build_backbone_batch(phi, psi)["CA"]
        got = pj.build_ca(phi, psi)[0]
        worst[n] = float(np.abs(got - ref).max())
    #: 1e-12 A is four orders below the 1e-8 A at which any downstream stage could see it,
    #: and the emitted coordinates do not come from here in any case.
    assert max(worst.values()) < 1e-12, worst
    assert worst[80] > worst[3], "error must accumulate with chain length, or it is a bug"


def test_scan_builder_reproduces_N_and_C_not_only_CA():
    """The gradient reads the frame axes, which are the N and C the objective never uses.
    If only CA agreed, the axes could be wrong and the gradient silently biased."""
    geo = core.backend("geometry")
    rng = np.random.default_rng(12)
    n = 17
    phi = rng.uniform(-math.pi, math.pi, (4, n))
    psi = rng.uniform(-math.pi, math.pi, (4, n))
    ref = geo.build_backbone_batch(phi, psi)
    G = pj.frames(phi, psi)
    assert np.abs(G[:, 1::3, :3, 3] - ref["N"][:, 1:]).max() < 1e-12
    assert np.abs(G[:, 2::3, :3, 3] - ref["CA"][:, 1:]).max() < 1e-12
    assert np.abs(G[:, 3::3, :3, 3] - ref["C"][:, 1:]).max() < 1e-12
    assert np.abs(G[:, 0, :3, 3] - ref["C"][:, 0]).max() < 1e-12


def test_hand_written_cross_is_bit_identical_to_numpy():
    """`numpy.cross` was 60% of the reference builder and 31% of the optimised gradient,
    almost all of it argument handling.  Replacing it is only legitimate if it is EXACTLY
    the same arithmetic: 0.0, not 'close'."""
    rng = np.random.default_rng(31)
    for shp in ((3,), (5, 3), (13, 3), (40, 3), (7, 11, 3)):
        a = rng.normal(size=shp)
        b = rng.normal(size=shp)
        assert np.abs(pj._cross3(a, b) - np.cross(a, b)).max() == 0.0, shp


def test_hand_written_det_agrees_in_sign_with_lapack():
    """It is used only for the reflection test, on a product of orthogonal matrices whose
    determinant is exactly +-1 -- nowhere near where a cofactor expansion and an LU could
    disagree.  Checked on rotations AND on reflections, since both occur."""
    rng = np.random.default_rng(32)
    for _ in range(200):
        Q = np.linalg.qr(rng.normal(size=(3, 3)))[0]
        for M in (Q, Q * np.array([-1.0, 1.0, 1.0])):
            assert abs(pj._det3(M) - np.linalg.det(M)) < 1e-13
            assert np.sign(pj._det3(M)) == np.sign(np.linalg.det(M))


def test_frames_are_rotations_so_the_rotation_axes_are_unit():
    """The gradient uses column 0 of every frame as a rotation axis and takes its length
    for granted.  If the scan drifted off SO(3) the gradient would be scaled wrong."""
    rng = np.random.default_rng(13)
    G = pj.frames(rng.uniform(-math.pi, math.pi, (3, 30)),
                  rng.uniform(-math.pi, math.pi, (3, 30)))
    R = G[..., :3, :3]
    I = np.einsum("...ji,...jk->...ik", R, R)
    assert np.abs(I - np.eye(3)).max() < 1e-13
    assert np.abs(np.linalg.det(R) - 1.0).max() < 1e-13


# =====================================================================================
# 2. THE ANALYTIC GRADIENT
# =====================================================================================
def test_analytic_gradient_matches_central_differences():
    """The mandate's check.  Central differences at h=1e-6 are themselves only good to
    ~1e-9 (truncation O(h^2 f''') against roundoff O(eps/h)), so 1e-8 is the floor of the
    comparison, not of the gradient."""
    rng = np.random.default_rng(7)
    worst = 0.0
    for n in (4, 9, 14, 25):
        C = rng.normal(0.0, 5.0, (n, 3))
        seq = "".join(rng.choice(list(pj.AA20), n))
        for lam, pen in ((0.0, None), (0.3, pj.RamaHingePenalty(seq, 0)),
                         (0.3, pj.RamaPenalty(seq, 0)), (0.3, pj.PhiPosPenalty(seq, 0))):
            x = rng.uniform(-math.pi, math.pi, 2 * n)
            fg = pj._make_fg(C, pen, lam, None, n, "analytic")
            ref_c = C - C.mean(0)
            ref_sq = float((ref_c ** 2).sum())

            def f(v):
                CA = pj.build_ca(v[None, :n], v[None, n:])[0][0]
                r = pj._rmsd_and_dca(CA, ref_c, ref_sq)[0]
                return r + (lam * float(pen(v[:n], v[n:])[0]) if lam and pen else 0.0)

            ga = fg(x)[1]
            h = 1e-6
            E = np.eye(2 * n) * h
            gn = np.array([(f(x + E[k]) - f(x - E[k])) / (2 * h) for k in range(2 * n)])
            worst = max(worst, float(np.abs(ga - gn).max()))
    assert worst < 1e-7, worst


def test_geometrically_inert_parameters_have_exactly_zero_gradient():
    """`phi[0]` and `psi[n-1]` are never read by the builder and are excluded from the
    penalty mask, so they are free parameters that cannot move the objective.  The
    reference's one-sided difference returns exactly 0.0 for them; so must this."""
    rng = np.random.default_rng(8)
    n = 12
    C = rng.normal(0.0, 5.0, (n, 3))
    seq = "".join(rng.choice(list(pj.AA20), n))
    pen = pj.RamaHingePenalty(seq, 2)
    x = rng.uniform(-math.pi, math.pi, 2 * n)
    for grad in pj.GRAD_MODES:
        g = pj._make_fg(C, pen, 0.3, None, n, grad)(x)[1]
        assert g[0] == 0.0, (grad, g[0])
        assert g[2 * n - 1] == 0.0, (grad, g[2 * n - 1])


def test_rmsd_matches_the_canonical_kabsch():
    """The inner loop scores with its own residual so it can also produce the rotation.
    It must be the SAME number the pipeline's canonical routine returns, or the objective
    being minimised is not the objective being reported."""
    aud = core.backend("numerics")
    rng = np.random.default_rng(9)
    for n in (5, 13, 30):
        C = rng.normal(0.0, 4.0, (n, 3))
        ref_c = C - C.mean(0)
        ref_sq = float((ref_c ** 2).sum())
        P = rng.normal(0.0, 4.0, (11, n, 3))
        want = aud.kabsch_rmsd_batch(P, C)
        assert np.abs(pj._rmsd_batch(P, ref_c, ref_sq) - want).max() < 5e-15
        got = np.array([pj._rmsd_and_dca(p, ref_c, ref_sq)[0] for p in P])
        assert np.abs(got - want).max() < 5e-15


def test_weighted_kabsch_reduces_to_the_unweighted_one():
    """Uniform weights must give the unweighted routine's residual exactly, or a weighted
    number could never be compared with the unweighted one it is supposed to bound."""
    aud = core.backend("numerics")
    rng = np.random.default_rng(10)
    n = 14
    C = rng.normal(0.0, 4.0, (n, 3))
    P = rng.normal(0.0, 4.0, (6, n, 3))
    w = np.ones(n)
    #: the weights are normalised to sum 1 inside, so the weighted residual is ALREADY a
    #: per-atom mean and the two routines are directly comparable with no 1/sqrt(n)
    assert np.abs(pj.wkabsch_rmsd_batch(P, C, w)
                  - aud.kabsch_rmsd_batch(P, C)).max() < 1e-14


# =====================================================================================
# 3. THE PENALTY IS THE SAME PENALTY
# =====================================================================================
def test_penalty_is_bit_identical_to_the_reference(inputs):
    """The tables are the trained prior.  Consolidating the code that reads them must not
    change a single value: exactly 0.0, on real sequences and real folds."""
    ref = pytest.importorskip("s8.project")
    rng = np.random.default_rng(4)
    worst = 0.0
    for row in inputs["rows"][:12]:
        for kind in ("ramah", "rama", "rama20", "phip"):
            a = pj.make_penalty(kind, row["seq"], int(row["fold"]))
            b = ref.make_penalty(kind, row["seq"], int(row["fold"]))
            n = row["n"]
            phi = rng.uniform(-math.pi, math.pi, (5, n))
            psi = rng.uniform(-math.pi, math.pi, (5, n))
            worst = max(worst, float(np.abs(a(phi, psi) - b(phi, psi)).max()))
            assert (a.mask == b.mask).all(), kind
    assert worst == 0.0, worst


def test_hinge_thresholds_are_bit_identical(inputs):
    ref = pytest.importorskip("s8.project")
    assert np.abs(pj.hinge_thresholds() - ref.hinge_thresholds()).max() == 0.0
    assert np.abs(pj.logp_tables("rama") - ref.logp_tables("rama")).max() == 0.0
    assert np.abs(pj.logp_tables("rama20") - ref.logp_tables("rama20")).max() == 0.0


# =====================================================================================
# 4. EQUIVALENCE ON REAL PIPELINE INPUTS
#
#    The bar here was RAISED after an independent audit.  It used to be "within optimiser
#    path noise"; it is now EXACTLY 0.0 for the shipped mode, because the audit showed that
#    a forward map agreeing to 1e-13 A still emits a different structure on 126/126 real
#    targets.  The projection is degenerate and L-BFGS-B starts far from any minimum, so
#    "numerically negligible" and "returns the same structure" are not the same claim.
# =====================================================================================
def test_the_shipped_objective_is_bit_identical_to_the_reference(inputs):
    """THE PROOF HALF.  L-BFGS-B is deterministic given (f, g), so if the shipped mode's
    objective and gradient are bit-identical to the reference's at arbitrary points, the
    trajectory and therefore the emitted structure must be identical too.

    The reference side is TRANSCRIBED from `s8.project.fit_prior` rather than called
    through it, so that what is being compared against is the reference's own arithmetic
    and not a wrapper that might already share code with the thing under test."""
    ref = pytest.importorskip("s8.project")
    audit = pytest.importorskip("s7.audit")
    geo = core.backend("geometry")
    rng = np.random.default_rng(99)
    wf = wg = 0.0
    for row in inputs["rows"][:8]:
        C = np.asarray(row["C"], float)
        n = len(C)
        pa = pj.make_penalty("ramah", row["seq"], int(row["fold"]))
        pb = ref.make_penalty("ramah", row["seq"], int(row["fold"]))
        for lam, qa, qb in ((0.0, None, None), (0.3, pa, pb)):
            E = np.eye(2 * n) * 1e-5

            def ref_fg(x, lam=lam, qb=qb, E=E, n=n, C=C):
                X = np.vstack([x[None], x[None] + E])
                B = geo.build_backbone_batch(X[:, :n], X[:, n:])["CA"]
                r = audit.kabsch_rmsd_batch(B, C)
                if lam and qb is not None:
                    r = r + lam * qb(X[:, :n], X[:, n:])
                return float(r[0]), (r[1:] - r[0]) / 1e-5

            fg = pj._make_fg(C, qa, lam, None, n, "exact")
            for _ in range(3):
                x = rng.uniform(-math.pi, math.pi, 2 * n)
                f1, g1 = ref_fg(x)
                f2, g2 = fg(x)
                wf = max(wf, abs(f1 - f2))
                wg = max(wg, float(np.abs(g1 - g2).max()))
    assert wf == TOL_EXACT, wf
    assert wg == TOL_EXACT, wg


def test_exact_builder_is_bit_identical_to_the_reference_builder():
    """And the forward map itself, at every batch width the objective actually uses."""
    geo = core.backend("geometry")
    rng = np.random.default_rng(77)
    for n in (2, 3, 9, 13, 16, 21, 40, 80):
        for B in (1, 3, 2 * n + 1):
            phi = rng.uniform(-math.pi, math.pi, (B, n))
            psi = rng.uniform(-math.pi, math.pi, (B, n))
            ref = geo.build_backbone_batch(phi, psi)["CA"]
            assert np.abs(pj.build_ca_exact(phi, psi) - ref).max() == TOL_EXACT, (n, B)


def test_live_subset_is_bit_identical_end_to_end(inputs):
    """THE MEASUREMENT HALF, live, so a stale `verify/project_exactness.json` cannot carry
    the suite: run `s8.project` and the shipped mode on real coordinate averages and
    compare the emitted coordinates and torsions themselves at exactly 0.0."""
    ref = pytest.importorskip("s8.project")
    for row in inputs["rows"][:3]:
        want = _run(ref, row)
        got = _run(pj, row, "exact")
        for k in (0, 1):                                  # the lam=0 fit and the lam=0.3 arm
            for j in (0, 1, 2):                           # CA, phi, psi
                assert np.abs(np.asarray(got[k][j]) - np.asarray(want[k][j])).max() \
                    == TOL_EXACT, (row["pdb"], k, j)
            assert got[k][3] == want[k][3], (row["pdb"], k, got[k][3], want[k][3])


def test_measured_equivalence_of_the_shipped_arm_over_all_126_targets(measured):
    """The four-arm table's own verdict on the shipped mode: emitted coordinates, torsions
    and objective identical to `s8.project` on every target, at exactly 0.0.

    IDENTITY IS READ OFF `darm_maxabs`, NOT off an RMSD.  A Kabsch RMSD between a structure
    and an identical copy of it is not 0.0 but ~1e-7: the residual is a cancellation of two
    large sums and the square root magnifies what survives.  Using it as the identity test
    would have hidden a real 1e-7 A difference and, worse, would have looked like a pass."""
    agg = measured["agg"]
    assert "ex" in agg, "regenerate with `python -m core.project equiv`"
    a = agg["ex"]
    assert a["darm_maxabs_worst"] == TOL_EXACT, a["darm_maxabs_worst"]
    assert a["dfval_worst"] == TOL_EXACT, a["dfval_worst"]
    assert a["drmsd_absmax"] == TOL_EXACT, a["drmsd_absmax"]
    assert a["n_moved_1e6"] == 0, a["n_moved_1e6"]
    #: and the aggregate the instrument is read on is the reference's, digit for digit
    assert a["synthesis"] == agg["ref"]["synthesis"], (a["synthesis"],
                                                       agg["ref"]["synthesis"])
    assert a["fit"] == agg["ref"]["fit"], (a["fit"], agg["ref"]["fit"])


def test_measured_exactness_over_all_126_targets():
    """The persisted whole-instrument check: eleven scalars per target, every one of them a
    function of the emitted coordinates, all required to be exactly 0.0 against the
    reference columns `s8.project` itself produced."""
    path = os.path.join(ROOT, "verify", "project_exactness.json")
    if not os.path.exists(path):
        pytest.skip("run `python -m core.project exactness`")
    with open(path) as fh:
        agg = json.load(fh)["agg"]
    assert agg["n_targets"] >= 120, agg["n_targets"]
    assert agg["not_identical"] == [], agg["not_identical"]
    assert agg["worst_abs_diff"] == TOL_EXACT, agg["worst_abs_diff"]
    assert agg["n_bit_identical"] == agg["n_targets"]


def test_the_fast_arm_is_priced_not_assumed(measured):
    """The `analytic` arm is a DIFFERENT PIPELINE and this test says so in numbers.

    It is not asserted equivalent -- it demonstrably is not -- but it may not silently move
    the aggregate the instrument is read on, and the per-target damage is pinned here so
    that a change to it shows up as a failure rather than as a quiet improvement."""
    agg = measured["agg"]
    assert agg["n_targets"] >= 120, agg["n_targets"]
    a = agg["an"]
    assert a["n_moved_1e6"] == agg["n_targets"], a["n_moved_1e6"]
    assert a["darm_rmsd_worst"] > 1.0, a["darm_rmsd_worst"]
    #: what it IS allowed to do: not move the aggregate, not be systematically worse
    assert abs(a["drmsd_mean"]) < TOL_MEAN_RMSD, a["drmsd_mean"]
    assert a["fval"] - agg["ref"]["fval"] < 5e-3, (a["fval"], agg["ref"]["fval"])


def test_the_reference_does_not_reproduce_itself_either():
    """THE CONTROL, and the reason the fast arm's differences are a property of the STAGE.

    A CA-RMSD is exactly invariant under a rigid motion of the target cloud and the penalty
    never sees the cloud, so `s8.project` run on C and on R.C + t is the same problem with
    every floating-point operation in the Kabsch differently rounded.  Where it returns a
    different structure, no reimplementation could have been "equivalent" except bit-exactly
    -- which is why the shipped mode is bit-exact and the fast one is keyed separately."""
    import glob
    #: the control is RAM-gated and stopped at 77 of 126 on this box, so it lives on the
    #: partial filename by design -- a partial table may not occupy the canonical path.
    #: 77 targets is far more than this assertion needs; the point is a property of the
    #: stage, not an aggregate over the whole instrument.
    cands = ([os.path.join(ROOT, "verify", "project_stability.json")]
             + sorted(glob.glob(os.path.join(ROOT, "verify",
                                             "project_stability_partial*.json"))))
    path = next((c for c in cands if os.path.exists(c)), None)
    if path is None:
        pytest.skip("run `python -m core.project stability`")
    with open(path) as fh:
        agg = json.load(fh)["agg"]
    assert agg["n_targets"] >= 20, agg["n_targets"]
    #: this must keep being TRUE, not be fixed: it is a measurement of the stage
    assert agg["darm_rmsd_worst"] > 1e-3, agg["darm_rmsd_worst"]


def test_the_instrument_numbers_are_reproduced(measured):
    """The numbers this stage is pinned to on the 126-target instrument.

    `EXPECT_FIT` is the multi-start lam=0 projection -- `s9/synth.py`'s own EXPECT and
    benchmark60's pre-registered `tuning_mean`.  The 3.2005 that circulated as the target
    is `s9.synth.fit_w` SINGLE-start from the medoid's torsions, a different construction,
    and pinning that one here would have made a correct stage look broken."""
    agg = measured["agg"]
    assert abs(agg["avg_rmsd"] - EXPECT_RMSD_AVG) < 1e-3, agg["avg_rmsd"]
    assert abs(agg["ref"]["fit"] - EXPECT_FIT) < 1e-3, agg["ref"]["fit"]
    assert abs(agg["ref"]["synthesis"] - EXPECT_SYNTHESIS) < 0.01, agg["ref"]["synthesis"]
    assert abs(agg["ref"]["projection_cost"] - EXPECT_PROJECTION_COST) < 0.01, \
        agg["ref"]["projection_cost"]


def test_the_optimisation_is_actually_faster(measured):
    """The whole point.  Measured on the same 126 real inputs in the same process as the
    reference it is divided by, so machine contention cancels out of the ratio."""
    agg = measured["agg"]
    assert agg["speedup_analytic"] > 5.0, agg["speedup_analytic"]
    assert agg["speedup_fd"] > 2.0, agg["speedup_fd"]
    path = os.path.join(ROOT, "verify", "project_exactness.json")
    if os.path.exists(path):
        with open(path) as fh:
            ex = json.load(fh)["agg"]
        #: and the SHIPPED, bit-exact mode is faster too -- that is what actually ships
        assert agg["s_per_target"]["ref"] / ex["s_per_target"] > 1.5, ex["s_per_target"]


# =====================================================================================
# 5. GEOMETRY VALIDITY ON EVERY EMISSION
# =====================================================================================
def test_every_emission_is_legal_geometry(measured):
    """Bond, pseudo-angle, clashes on every emission of every arm.

    TWO DIFFERENT BARS, and conflating them was a real mistake in the first version of this
    file.  The SHIPPED arm is held to EQUIVALENCE -- every geometry statistic identical to
    the reference's, because it is the same structures.  The fast arms are different
    pipelines, so they are held to LEGALITY (bond, bond SD, no clash epidemic, L-handed)
    and their differences are RECORDED rather than forbidden.

    Asserting `clashes <= ref` on the fast arms was the wrong shape and it failed
    honestly: `fd` returns a sub-4 A non-local contact on 3 targets where the reference has
    2.  That is a fact about landing on a different branch, not a defect to tolerate away,
    and on an arm already known to move 126/126 structures it is not illegality either.  A
    CA-RMSD gain on an illegal chain is still not a gain -- which is why the absolute bond
    and chirality bars below apply to every arm without exception.
    """
    agg = measured["agg"]
    for tag in ARMS:
        a = agg[tag]
        #: absolute legality, every arm: the builder places every bond from a constant
        assert abs(a["step_mean"] - EXPECT_BOND) < 5e-3, (tag, a["step_mean"])
        assert a["step_sd_max"] < 1e-14, (tag, a["step_sd_max"])
        assert a["min_nonlocal_CA"] > 2.5, (tag, a["min_nonlocal_CA"])
        #: the projection returns a somewhat over-helical chain against the natives'
        #: 103.99 -- pinned as a FACT about the stage, not as a target
        assert 88.0 < a["pseudoangle_mean"] < NATIVE_PSEUDOANGLE, \
            (tag, a["pseudoangle_mean"])

    #: EQUIVALENCE, shipped arm only: the whole table, digit for digit
    for k in ("step_mean", "step_sd_max", "pseudoangle_mean", "posphi_con_nongly",
              "min_nonlocal_CA", "clashes", "L_worst", "synthesis", "fit"):
        assert agg["ex"][k] == agg["ref"][k], (k, agg["ex"][k], agg["ref"][k])

    #: PRICED, fast arms: recorded so a change fails a test rather than passing quietly
    assert agg["ref"]["clashes"] == 2.0, agg["ref"]["clashes"]
    assert agg["fd"]["clashes"] == 3.0, agg["fd"]["clashes"]
    assert agg["an"]["clashes"] == 3.0, agg["an"]["clashes"]


def test_positive_phi_on_the_constrained_non_glycine_set(measured):
    """`phi[0]` is never read by the builder, so including it measures a coin flip; the
    honest number is over `phi[1:]`, non-glycine.  The `ramah` prior puts it in the
    2.90-5.66% band real backbones occupy, and the optimised arms must keep it there."""
    agg = measured["agg"]
    for tag in ("ref", "fd", "an"):
        r = agg[tag]["posphi_con_nongly"]
        assert r < 0.10, (tag, r)
        assert abs(r - agg["ref"]["posphi_con_nongly"]) < 0.01, (tag, r)


def test_every_emission_is_L_handed(measured):
    """THE CHIRALITY ASSERTION.  A distance objective is exactly mirror-blind and its
    lowest-objective multi-start selects enantiomers; this objective is a COORDINATE
    distance with reflections forbidden in the superposition, so the branch it picks is
    chirality-resolved.  That safety is a property of the objective and would be silently
    lost if the objective were ever reformulated, so it is asserted on the emission."""
    agg = measured["agg"]
    for tag in ("ref", "fd", "an"):
        assert agg[tag]["L_worst"] == 1.0, (tag, agg[tag]["L_worst"])


def test_L_signature_detects_a_mirror():
    """A test that only ever returns +1 is not a test.  Reflect the chain and it must flip,
    or the assertion above is measuring nothing."""
    geo = core.backend("geometry")
    rng = np.random.default_rng(21)
    n = 15
    phi = rng.uniform(-math.pi, math.pi, n)
    psi = rng.uniform(-math.pi, math.pi, n)
    assert pj.l_signature(phi, psi) == 1.0
    B = geo.build_backbone_batch(phi[None], psi[None])
    N, CA, C, CB = (B[k][0] * np.array([-1.0, 1.0, 1.0]) for k in ("N", "CA", "C", "CB"))
    v = ((N - CA) * np.cross(C - CA, CB - CA)).sum(1)
    assert float(np.sign(v).mean()) == -1.0


# =====================================================================================
# 6. THE PARAMETERS THAT ARE SCIENCE, NOT TUNING
# =====================================================================================
def test_the_starts_are_the_same_four():
    """Multi-start is load-bearing: the projection is degenerate, warm-starting cannot
    cross between the two branches, and reducing the starts would be weakening the science
    rather than optimising it.  The four are `s8.consensus2.FIT_STARTS`, by value."""
    cc = pytest.importorskip("s8.consensus2")
    assert tuple(pj.STARTS) == tuple(cc.FIT_STARTS)
    assert len(pj.STARTS) == 4


def test_the_pinned_constants_are_unchanged():
    ref = pytest.importorskip("s8.project")
    geo = core.backend("geometry")
    assert pj.FD_EPS == 1e-5
    assert (pj.RB, pj.SIGMA_BINS, pj.PSEUDO) == (ref.RB, ref.SIGMA_BINS, ref.PSEUDO)
    assert pj.HINGE_PCT == ref.HINGE_PCT and pj.PHIP_TAU == ref.PHIP_TAU
    assert pj.CLASSES == ref.CLASSES and pj.AA20 == ref.AA20
    for k in ("BOND_N_CA", "BOND_CA_C", "BOND_C_N", "ANGLE_N_CA_C", "ANGLE_CA_C_N",
              "ANGLE_C_N_CA", "OMEGA_TRANS"):
        assert getattr(pj, k) == getattr(geo, k), k


def test_lam_path_runs_the_declared_number_of_solves(inputs):
    """`multi=True` at production settings is 4 generic starts at lam=0, one continuation
    at lam=0.3, and 4 generic starts again at lam=0.3 -- nine L-BFGS-B solves.  A "speed
    optimisation" that quietly dropped one would be a change of method."""
    row = inputs["rows"][0]
    seen = []
    orig = pj.fit_prior

    def spy(*a, **k):
        seen.append(k.get("maxiter"))
        return orig(*a, **k)

    pj.fit_prior = spy
    try:
        pen = pj.make_penalty("ramah", row["seq"], int(row["fold"]))
        pj.lam_path(np.asarray(row["C"], float), pen, (0.0, 0.3), maxiter=300, multi=True)
    finally:
        pj.fit_prior = orig
    assert len(seen) == 9, seen
    assert set(seen) == {300}, seen


def test_the_backend_switch_routes_here():
    """A consolidated module going live is the act of committing the file; this asserts the
    contract is complete so the switch cannot half-fire."""
    assert core.backend_name("project") == "core.project"
    for sym in core.CONTRACT["project"]:
        assert hasattr(pj, sym), sym
