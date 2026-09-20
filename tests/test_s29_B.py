"""Tests for `s29/s29_B_compat.py` (S29 lane B: the compatibility Hamiltonian H = diag(E) - J M).

What each test pins:
  * the matrix builders: symmetry, zero diagonal where required, unit spectral norm, zero
    padding rows/columns;
  * the CENTERING IDENTITIES that the whole lane rests on -- A_c 1 = 0 and G 1 = 0 -- and their
    corollary, the SIGN-MIXING LEMMA: every eigenvector of a centered matrix with lambda != 0 has
    entries of both signs (prereg section 3.2);
  * G's construction: G = Delta Delta^T / n_res is PSD, its rank is at most 3 n_res - 3, and its
    nonzero spectrum equals that of the pool's deviation covariance Delta^T Delta / n_res
    (the same nonzero spectrum, which is what "G is the operator whose eigenvectors are the
    pool's shape modes" means);
  * the stable-rank bound r_stable(G) <= rank(G) <= 3 n_res - 3 <= 45, which is why lane T's law
    (Var = r_stable / D^2) caps this family's gradient variance by construction;
  * the two controls: PERM and SPEC carry M's spectrum to 1e-10, PERM keeps the entry multiset,
    SPEC destroys the eigenvector structure;
  * the INDEPENDENT gradient estimator: parameter shift equals central finite differences, and
    equals `s27.s28_B_hop.grad_hop_paramshift` (so measurement 1's independence is in the
    matrices and the driver, and any disagreement is a real bug, not a convention);
  * the exact ground state: one-hot argmin at J = 0, the eigen-equation at J != 0, and the
    global-sign convention being well defined for a centered M (where `sum(v) ~ 0`);
  * the signed readout R4 reduces to the uniform average on a flat state over a set;
  * NaN-poison: the native-free half is bit-identical and every ORACLE column is NaN.

Run:  python -m pytest tests/test_s29_B.py -q -p no:cacheprovider
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from core import quantum as Q                       # noqa: E402
from s29 import s29_B_compat as C                   # noqa: E402


def _rand_pool(k=40, n_res=9, seed=5):
    rng = np.random.default_rng(seed)
    return rng.normal(0.0, 3.0, (k, n_res, 3)), rng


def _pairwise(W):
    import core
    kb = core.backend("numerics").kabsch_rmsd_batch
    W = np.asarray(W, float)
    D = np.empty((len(W), len(W)))
    for a in range(len(W)):
        D[a] = kb(W, W[a])
    D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    return D


# ------------------------------------------------------------------------ the matrix builders
def test_gaussian_similarity_reproduces_the_s28_graph():
    """The independent builder is the S28 graph, or measurement 1 is not comparable to S28."""
    from s27 import s28_B_hop as B
    W, _ = _rand_pool()
    D = _pairwise(W)
    mine = C.gaussian_similarity(D)
    theirs = B.kernel_graph(D)
    assert abs(mine["sigma"] - theirs["sigma"]) < 1e-12
    assert np.max(np.abs(mine["A"] - theirs["A"])) < 1e-12


def test_every_matrix_is_symmetric_and_unit_spectral_norm():
    W, rng = _rand_pool()
    D = _pairwise(W)
    A = C.gaussian_similarity(D)["A"]
    Ac = C.unit_spectral(C.double_center(A))
    G = C.unit_spectral(C.agreement_gram(W))
    for M in (A, Ac, G):
        assert np.max(np.abs(M - M.T)) < 1e-12
        assert abs(float(np.abs(np.linalg.eigvalsh(M)).max()) - 1.0) < 1e-10
    assert np.max(np.abs(np.diag(A))) == 0.0            # the similarity graph has no self-loop


def test_padding_rows_and_columns_are_zero():
    W, _ = _rand_pool(k=40)
    A = C.gaussian_similarity(_pairwise(W))["A"]
    P = C.pad(A, 64)
    assert P.shape == (64, 64)
    assert np.max(np.abs(P[:40, :40] - A)) == 0.0
    assert np.max(np.abs(P[40:, :])) == 0.0 and np.max(np.abs(P[:, 40:])) == 0.0


# ----------------------------------------------- the centering identities and the sign lemma
def test_centering_identities_A_c_one_is_zero_and_G_one_is_zero():
    W, _ = _rand_pool()
    A = C.gaussian_similarity(_pairwise(W))["A"]
    Ac = C.unit_spectral(C.double_center(A))
    G = C.unit_spectral(C.agreement_gram(W))
    one = np.ones(len(A))
    assert np.max(np.abs(Ac @ one)) < 1e-10
    assert np.max(np.abs(G @ one)) < 1e-10
    # ... and A itself does NOT annihilate 1 (it is the typicality projector S28 measured)
    assert np.max(np.abs(A @ one)) > 0.1


def test_sign_mixing_lemma_every_nonzero_eigenvector_has_both_signs():
    """The lemma the whole lane rests on: M 1 = 0 => every eigenvector with lambda != 0 is
    orthogonal to 1 => it has entries of both signs => a p-readout averages the two poles."""
    W, _ = _rand_pool()
    A = C.gaussian_similarity(_pairwise(W))["A"]
    for M in (C.unit_spectral(C.double_center(A)), C.unit_spectral(C.agreement_gram(W))):
        w, V = np.linalg.eigh(M)
        for j in range(len(w)):
            if abs(w[j]) > 1e-8:
                v = V[:, j]
                assert abs(float(v.sum())) < 1e-8 * len(v)
                assert (v > 1e-12).any() and (v < -1e-12).any()


def test_G_is_psd_low_rank_and_shares_the_pool_covariance_spectrum():
    """The rank cap is what makes lane T's law (Var = r_stable/D^2) bind this family.

    The EXACT bound on POSED windows is 3 n_res - 3: Kabsch centres every member on the
    reference's centroid, so the three translation directions are removed identically, while the
    three rotations are removed only to first order (hence the 3 n_res - 6 quoted in S29-L11 and
    in this lane's prereg, which is the practical figure, not an identity). Either way the cap is
    at most 45 on this instrument against the 266 a within-30x gradient variance would need.
    """
    from s12 import instrument as I
    for n_res in (9, 12, 16):
        W0, _ = _rand_pool(k=60, n_res=n_res, seed=n_res)
        W = I.superpose_batch(W0, W0[0])                         # the frame the readout averages in
        Delta = C.deviations(W)
        G = C.agreement_gram(W)
        w = np.linalg.eigvalsh(G)
        assert w.min() > -1e-9                                   # PSD
        rank = int((w > 1e-8 * w.max()).sum())
        assert rank <= 3 * n_res - 3                             # the exact cap (translation)
        assert rank <= 45                                        # the instrument's own ceiling
        # the same nonzero spectrum as the pool's deviation covariance
        S = (Delta.T @ Delta) / float(n_res)
        ws = np.sort(np.linalg.eigvalsh(S))[::-1][:rank]
        wg = np.sort(w)[::-1][:rank]
        assert np.max(np.abs(ws - wg)) < 1e-8 * max(1.0, wg.max())
        # and therefore the stable rank is capped, which is lane T's law applied to this family
        Gn = C.unit_spectral(G)
        assert float((Gn ** 2).sum()) <= rank + 1e-6


# ------------------------------------------------------------------------------ the controls
def test_perm_and_spec_controls_carry_the_spectrum():
    W, rng = _rand_pool()
    A = C.gaussian_similarity(_pairwise(W))["A"]
    Ac = C.unit_spectral(C.double_center(A))
    w0 = np.sort(np.linalg.eigvalsh(Ac))
    P = C.perm_matrix(Ac, rng.permutation(len(Ac)))
    S = C.spec_matrix(Ac, np.random.default_rng(11))
    for M in (P, S):
        assert np.max(np.abs(M - M.T)) < 1e-10
        assert np.max(np.abs(np.sort(np.linalg.eigvalsh(M)) - w0)) < 1e-10
    # PERM keeps the entry multiset; SPEC does not (it keeps only the spectrum)
    assert np.allclose(np.sort(P.reshape(-1)), np.sort(Ac.reshape(-1)))
    assert not np.allclose(np.sort(S.reshape(-1)), np.sort(Ac.reshape(-1)))
    # SPEC still annihilates nothing in particular: its kernel is rotated away from 1
    assert np.max(np.abs(S @ np.ones(len(S)))) > 1e-8


# -------------------------------------------------------------- the independent hop gradient
def test_hop_gradient_paramshift_equals_finite_differences():
    rng = np.random.default_rng(3)
    circ = Q.StatevectorCircuit(5, 3)
    M = rng.random((circ.dim, circ.dim))
    M = 0.5 * (M + M.T)
    for _ in range(3):
        th = rng.normal(0.0, 0.6, circ.n_params())
        g_ps = C.hop_grad_paramshift(circ, th, M)
        g_fd = C.hop_grad_fd(circ, th, M)
        assert g_ps.shape == (circ.n_params(),)
        assert np.max(np.abs(g_ps - g_fd)) < 1e-6


def test_hop_gradient_agrees_with_the_s28_implementation():
    """Measurement 1's independence is in the MATRICES and the driver; the shift rule itself
    must agree with S28's or a disagreement in the rows would be a convention, not a bug."""
    from s27 import s28_B_hop as B
    rng = np.random.default_rng(7)
    circ = Q.StatevectorCircuit(4, 3)
    M = rng.random((circ.dim, circ.dim))
    M = 0.5 * (M + M.T)
    th = rng.normal(0.0, 0.6, circ.n_params())
    assert np.max(np.abs(C.hop_grad_paramshift(circ, th, M)
                         - B.grad_hop_paramshift(circ, th, M))) < 1e-12


def test_hop_variance_returns_the_stable_rank_and_its_prediction():
    W, _ = _rand_pool(k=16)
    A = C.gaussian_similarity(_pairwise(W))["A"]
    r = C.hop_variance(4, A, n_theta=12, fd_check=True)
    assert r["dim"] == 16 and r["P"] == 12
    assert abs(r["r_stable"] - float((A ** 2).sum())) < 1e-9      # unit spectral norm
    assert abs(r["pred_var"] - r["fro2"] / 16 ** 2) < 1e-15
    assert r["fd_rel_max"] < 1e-5
    assert r["var_g0"] > 0


# --------------------------------------------------------------------------- the ground state
def test_ground_state_is_the_argmin_at_J_zero_and_solves_the_eigen_equation():
    rng = np.random.default_rng(13)
    dim = 32
    E = rng.normal(0.0, 1.0, dim)
    M = rng.random((dim, dim))
    M = 0.5 * (M + M.T)
    np.fill_diagonal(M, 0.0)
    gs0 = C.ground_state(E, M, 0.0)
    assert int(np.argmax(gs0["p"])) == int(np.argmin(E))
    assert gs0["p"].max() > 1 - 1e-12
    for J in (-1.0, 1.0):
        gs = C.ground_state(E, M, J)
        Hm = np.diag(E) - J * M
        assert np.max(np.abs(Hm @ gs["psi"] - gs["e0"] * gs["psi"])) < 1e-9
        assert abs(float(gs["p"].sum()) - 1.0) < 1e-12


def test_ground_state_sign_convention_is_well_defined_for_a_centered_matrix():
    """S28 fixed the global sign by `sum(v) >= 0`; for a centered M that sum is ~0 and its sign
    is noise. The largest-|entry| convention must be deterministic under a sign flip of eigh."""
    W, _ = _rand_pool(k=32)
    A = C.gaussian_similarity(_pairwise(W))["A"]
    Ac = C.unit_spectral(C.double_center(A))
    E = np.zeros(len(Ac))
    gs = C.ground_state(E, Ac, 5.0)
    v = gs["psi"]
    assert v[int(np.argmax(np.abs(v)))] > 0
    assert abs(float(v.sum())) < 1e-6 * np.abs(v).sum() + 1e-9   # the pole symmetry, measured


# ------------------------------------------------------------------------ the signed readout
def test_signed_readout_reduces_to_the_uniform_average_on_a_flat_state():
    from s27 import s28_A_amp as AMP
    W, rng = _rand_pool(k=40)
    top = np.sort(rng.choice(40, 20, replace=False))
    frame = AMP.Frame(W, top)
    psi = np.zeros(40)
    psi[top] = 1.0
    Cmix, denom, w = AMP.readout(psi, frame)
    assert abs(denom - 20.0) < 1e-12
    assert np.max(np.abs(Cmix - frame.Wf[top].mean(0).reshape(frame.n, 3))) < 1e-10


# ---------------------------------------------------------------------------- NaN-poison
def test_selection_is_native_free_and_nan_poison_is_bit_identical(monkeypatch):
    from s24 import d_harness as H
    from s27 import run_pool as RP
    from s29 import s29_B_compat as CC
    rng = np.random.default_rng(23)
    k, n_res = 40, 9
    W = rng.normal(0.0, 3.0, (k, n_res, 3))
    E = (np.arange(k, dtype=float) - (k - 1) / 2.0)[rng.permutation(k)] / k
    key = RP.rng_for("SYN1", "tiekey").random(k)
    mats = ("A_c", "G")
    jg = (-1.0, 1.0)
    a = CC.gs_target_selection(W, E, key, "SYN1", j_grid=jg, mats=mats)
    b = CC.gs_target_selection(W, E, key, "SYN1", j_grid=jg, mats=mats)
    assert len(a["cells"]) == len(mats) * len(jg)
    for ca, cb in zip(a["cells"], b["cells"]):
        assert ca["R3"]["idx"] == cb["R3"]["idx"]
        for R in ("R2", "R3", "R4"):
            # assert_array_equal treats NaNs in the same positions as equal; R4 is legitimately
            # undefined for a centered M (the pole symmetry -- that IS the lemma, measured here)
            np.testing.assert_array_equal(np.asarray(ca[R]["C"]), np.asarray(cb[R]["C"]))
        assert np.isfinite(np.asarray(ca["R2"]["C"])).all()
        assert np.isfinite(np.asarray(ca["R3"]["C"])).all()
        assert (ca["eta"] == cb["eta"] or (np.isnan(ca["eta"]) and np.isnan(cb["eta"])))
        assert ca["pr"] == cb["pr"]
    nat = rng.normal(0.0, 3.0, (n_res, 3))
    ok = H.Candidates(pdb="SYN1", n=n_res, seq="A" * n_res, fold=0, W=W, nat_ca=nat,
                      oracle_rr=np.ones(k))
    bad = H.Candidates(pdb="SYN1", n=n_res, seq="A" * n_res, fold=0, W=W,
                       nat_ca=np.full((n_res, 3), np.nan), oracle_rr=np.full(k, np.nan))
    calls = {"i": 0}

    def fake_channels(pdb):
        c = ok if calls["i"] == 0 else bad
        calls["i"] += 1
        return c, {"DIS": E}, {}

    monkeypatch.setattr(RP, "channels_for", fake_channels)
    monkeypatch.setattr(RP, "zr", lambda x: np.asarray(x, float))
    monkeypatch.setattr(CC, "MAT_ARMS", mats)
    monkeypatch.setattr(CC, "J_GRID", jg)
    ra = CC.gs_target(pdb="SYN1")
    rb = CC.gs_target(pdb="SYN1")
    assert len(ra) == len(rb) == len(mats) * len(jg)
    for x, y in zip(ra, rb):
        assert x["R3"]["idx"] == y["R3"]["idx"]
        assert (x["eta"] == y["eta"]) or (np.isnan(x["eta"]) and np.isnan(y["eta"]))
        assert x["pr"] == y["pr"] and x["sign_coh"] == y["sign_coh"]
        for R in ("R2", "R3", "R4"):
            assert np.isnan(y[R]["rmsd"])
        assert np.isfinite(x["rmsd_prod"]) and np.isnan(y["rmsd_prod"])
        assert "sign_correct" not in y or np.isnan(y.get("rmsd_pc1_bestsign", np.nan))


def test_ties_never_break_by_array_order_in_the_R3_readout():
    """An exactly tied probability vector must not read the pool's own order (project memory)."""
    from s27 import run_pool as RP
    rng = np.random.default_rng(31)
    k = 40
    p = np.full(k, 1.0 / k)
    key_a = RP.rng_for("SYN_A", "tiekey").random(k)
    key_b = RP.rng_for("SYN_B", "tiekey").random(k)
    ia = np.sort(RP.topm(-p, 10, key_a))
    ib = np.sort(RP.topm(-p, 10, key_b))
    assert not np.array_equal(ia, ib)                 # different keys -> different tied choices
    assert not np.array_equal(ia, np.arange(10))      # and neither is the array order
    del rng


# ================================================== measurement 5: tail-then-aggregate (TTA)
def test_tail_lambda_is_the_deployed_cvar_tail_and_sums_to_one():
    from s29 import s29_B_tta as TT
    rng = np.random.default_rng(17)
    D = 64
    E = np.sort(rng.normal(0.0, 1.0, D))
    p = rng.random(D)
    p = p / p.sum()
    lam, strict, x_q, q, mass = TT.tail_lambda(E, p, 0.18)
    assert abs(float(lam.sum()) - 1.0) < 1e-12
    assert (lam >= 0).all()
    # the tail's E-average equals the exact CVaR value
    v, q2, _ = Q.cvar_exact(E, p, 0.18)
    assert abs(float(lam @ E) - v) < 1e-10 and abs(q - q2) < 1e-12
    # the strict tail carries p/alpha and everything above the VaR carries nothing
    assert np.max(np.abs(lam[strict] - p[strict] / 0.18)) < 1e-12
    above = np.array([y for y in range(D) if y not in set(strict.tolist()) and y != x_q])
    assert np.max(np.abs(lam[above])) == 0.0


def test_tta_gradient_matches_finite_differences():
    """The envelope gradient dR/dp = (W_y - W_xq)/alpha, chained through parameter shift."""
    from s27 import s28_A_amp as AMP
    from s29 import s29_B_tta as TT
    rng = np.random.default_rng(19)
    n_qubits, n_res = 5, 9
    circ = Q.StatevectorCircuit(n_qubits, 3)
    D = circ.dim
    E = np.sort(rng.normal(0.0, 1.0, D))
    Wf = rng.normal(0.0, 3.0, (D, 3 * n_res))

    class _Sur:                                   # a smooth surrogate with an exact gradient
        def value_grad(self, C):
            C = np.asarray(C, float)
            return float((C ** 2).sum()), 2.0 * C

    sur = _Sur()
    th = rng.normal(0.0, 0.6, circ.n_params())
    lam_f = 0.7
    val, g, p, _ = TT.objective_theta(circ, th, E, 0.18, 0.5, lam_f, Wf, sur, n_res)
    h = 1e-6
    gfd = np.zeros_like(g)
    for kk in range(len(th)):
        tp, tm = th.copy(), th.copy()
        tp[kk] += h
        tm[kk] -= h
        vp = TT.objective_theta(circ, tp, E, 0.18, 0.5, lam_f, Wf, sur, n_res)[0]
        vm = TT.objective_theta(circ, tm, E, 0.18, 0.5, lam_f, Wf, sur, n_res)[0]
        gfd[kk] = (vp - vm) / (2 * h)
    assert np.max(np.abs(g - gfd)) < 1e-4 * max(1.0, float(np.abs(g).max()))
    del AMP, val


def test_flat_report_reproduces_the_derived_flat_set():
    """Prereg B2.3: the f term is flat on EXACTLY the CVaR term's 437 directions, so the raw flat
    fraction is unchanged and the OVERLAP is what moves."""
    from s29 import s29_B_tta as TT
    rng = np.random.default_rng(23)
    n_res, D = 9, 64
    E = np.sort(rng.normal(0.0, 1.0, D))
    Wf = rng.normal(0.0, 3.0, (D, 3 * n_res))
    p = rng.random(D) ** 3
    p = p / p.sum()

    class _Sur:
        def value_grad(self, C):
            C = np.asarray(C, float)
            return float((C ** 2).sum()), 2.0 * C

    r0 = TT.flat_report(E, p, 0.18, 0.5, 0.0, Wf, _Sur(), n_res)
    r1 = TT.flat_report(E, p, 0.18, 0.5, 1.0, Wf, _Sur(), n_res)
    assert r0["flat_cvar"] == r1["flat_cvar"]
    assert r1["flat_f"] == r1["flat_cvar"]              # the same flat subspace, derived
    assert r1["flat_info_bearing"] == r1["flat_cvar"]
    assert r1["flat_readout_tta"] == r1["flat_cvar"]    # the readout's flat set is the same
    assert r1["overlap_tta"] == 0.0                     # every moving direction is seen
    assert r1["flat_readout_deployed"] == 1.0           # the deployed readout moves nowhere
    assert r1["nmove_readout_tta"] == r1["m_strict"]
    assert 0.0 < r1["flat_cvar"] < 1.0


def test_argmin_untied_never_reads_array_order_on_a_tie():
    """The project's named failure mode (`tie-breaking-leaks-the-pool-order`; lane D's S29-L38
    found the first version of the subset search taking `tie[0]`, which on a DIS-sorted pool is
    biased toward the very prefix the experiment tests)."""
    from s29 import s29_B_tta as TT
    v = np.zeros(50)                       # every index tied
    seen = set()
    for s_ in range(200):
        i, n = TT.argmin_untied(v, np.random.default_rng(s_))
        seen.add(i)
        assert n == 50
    assert len(seen) > 20                  # it spreads over the tie set, not index 0
    assert seen != {0}
    # a unique minimum is still returned exactly, with tie set size 1
    v2 = np.arange(50, dtype=float)
    v2[37] = -1.0
    i, n = TT.argmin_untied(v2, np.random.default_rng(0))
    assert i == 37 and n == 1
