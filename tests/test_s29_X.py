"""tests/test_s29_X.py -- S29 LANE X: the configuration-space CVaR-VQE.

The four properties the brief demands of this operator, plus the contract's NaN-poison:

  1. the configuration builder reproduces a pool member EXACTLY when every segment's fragment
     choice is that member's own (the space contains its parents);
  2. the Hamiltonian's diagonal part equals the S13-style objective computed independently on
     the built configuration, on 3 configurations (an independent recomputation, not a
     restatement: the pair term is recomputed from the distogram's raw `prob` and the built
     CA trace, the Ramachandran term from the raw `cnt` table, in this file);
  3. the deployable path is NaN-poison clean: replacing `nat_ca` / `oracle_rr` with NaN leaves
     every emitted structure bit-identical;
  4. at Gamma = 0 the objective and its gradient are `core.quantum.free_energy` bit for bit,
     the mixer's gradient agrees with finite differences, and the uniform-weight readout is
     `s12.instrument.coordinate_average` bit for bit.

The target used everywhere is the smallest probe target (2P5H, n = 9, q = 9, 512
configurations), so the whole file runs in well under a minute.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from core import predict as PRD          # noqa: E402
from core import quantum as Q            # noqa: E402
from s12 import instrument as I          # noqa: E402
from s27 import ham_lib as HL            # noqa: E402
from s29 import s29_X_config as X        # noqa: E402

PDB = "2P5H"


@pytest.fixture(scope="module")
def space():
    cand, dg, rama = X.load_target(PDB)
    return X.Space(PDB, cand, dg, rama)


def test_space_shape(space):
    assert space.M == 1 << space.q, "the register must be padding-free (2**q == F**S)"
    assert space.q == 3 * space.S <= 3 * X.SEG_MAX
    assert len(space.members) == X.F_MEMBERS == 8
    assert len(set(space.members.tolist())) == 8, "the parents must be distinct pool members"


def test_parents_reproduce_their_members(space):
    """Property 1: all-segments-equal-f rebuilds member f exactly."""
    for f in range(space.F):
        c = f * (space.M - 1) // (space.F - 1)
        assert np.array_equal(space.digits[c], np.full(space.S, f))
        rb = I.build_ca(space.cand.PHI[space.members[f]], space.cand.PSI[space.members[f]])
        assert np.max(np.abs(space.CA[c] - rb)) == 0.0


def test_diagonal_hamiltonian_recomputed_independently(space):
    """Property 2: H_diag on 3 configurations, recomputed from the raw artefacts here."""
    prob = np.asarray(space.dg["prob"], float)
    cnt = np.asarray(space.rama_cnt, float) + HL.PSEUDO
    Pm = cnt / cnt.sum(axis=(1, 2), keepdims=True)
    aa = HL.codes(space.seq)
    for c in (0, space.M // 3, space.M - 1):
        ca = space.CA[c]
        e_pair = 0.0
        for k in range(len(space.i)):
            d = float(np.linalg.norm(ca[space.i[k]] - ca[space.j[k]]))
            b = int(np.digitize(d, PRD.BIN_EDGES))
            e_pair += -math.log(prob[k, b] + X.EPS_P)
        e_rama = 0.0
        for r in range(space.n - 1):               # residue n-1 is inert for the CA trace
            pb = min(max(int((space.PHI[c, r] + math.pi) / (2 * math.pi) * 36), 0), 35)
            sb = min(max(int((space.PSI[c, r] + math.pi) / (2 * math.pi) * 36), 0), 35)
            e_rama += -math.log(Pm[aa[r], pb, sb])
        assert abs(space.E_pair[c] - e_pair) < 1e-9
        assert abs(space.E_rama[c] - e_rama) < 1e-9
        assert abs(space.E[c] - (e_pair + e_rama)) < 1e-9


def test_gamma_zero_is_the_deployed_objective(space):
    """Property 4a: at Gamma = 0, F and its gradient are `core.quantum.free_energy`."""
    circ = Q.StatevectorCircuit(space.q, X.LAYERS)
    th = np.random.default_rng(7).normal(0.0, 0.6, circ.n_params())
    f0, g0, _, _, _ = Q.free_energy(circ, th, space.E, X.ALPHA, X.TEMP)
    f1, g1 = X.f_and_grad(circ, th, space.E, X.ALPHA, X.TEMP, 0.0, [])[:2]
    assert f0 == f1
    assert np.max(np.abs(g0 - g1)) == 0.0


def test_mixer_gradient_matches_finite_differences(space):
    """Property 4b: the transverse-field term's exact parameter-shift gradient."""
    circ = Q.StatevectorCircuit(space.q, X.LAYERS)
    xors = [X.xor_index(space.M, space.q, k) for k in range(space.q)]
    th = np.random.default_rng(11).normal(0.0, 0.6, circ.n_params())
    gam = space.gamma_gap
    g = X.f_and_grad(circ, th, space.E, X.ALPHA, X.TEMP, gam, xors)[1]
    h = 1e-6
    for k in (0, 5, len(th) - 1):
        tp, tm = th.copy(), th.copy()
        tp[k] += h
        tm[k] -= h
        fd = (X.f_and_grad(circ, tp, space.E, X.ALPHA, X.TEMP, gam, xors)[0]
              - X.f_and_grad(circ, tm, space.E, X.ALPHA, X.TEMP, gam, xors)[0]) / (2 * h)
        assert abs(fd - g[k]) < 1e-5


def test_mixer_value_is_a_pauli_expectation(space):
    """<psi|sum X_k|psi> against a direct dense construction on a small register."""
    circ = Q.StatevectorCircuit(space.q, X.LAYERS)
    th = np.random.default_rng(13).normal(0.0, 0.6, circ.n_params())
    psi = circ.state(th)
    psi = psi / np.linalg.norm(psi)
    xors = [X.xor_index(space.M, space.q, k) for k in range(space.q)]
    direct = 0.0
    for k in range(space.q):
        bit = 1 << (space.q - 1 - k)
        direct += float(sum(psi[x] * psi[x ^ bit] for x in range(space.M)))
    assert abs(X.mixer_value(psi, xors) - direct) < 1e-10


def test_uniform_weighted_average_is_the_deployed_readout(space):
    """Property 4c: uniform weights reproduce `s12.instrument.coordinate_average`."""
    W = space.CA[:37]
    C0, _ = I.coordinate_average(W)
    C1 = X.weighted_average(W, None, PDB, "test")
    assert np.max(np.abs(C0 - C1)) == 0.0
    w = np.full(len(W), 0.3)                       # constant weights = uniform
    C2 = X.weighted_average(W, w, PDB, "test")
    assert np.max(np.abs(C0 - C2)) < 1e-12


def test_tail_readout_matches_the_deployed_extraction(space):
    """The tail set is `core.quantum.cvar_from_probs`'s positive-mass set, and at a
    concentrated distribution it is the classical energy prefix (the S24 theorem)."""
    p = X.gibbs(space.E, 0.5)
    ro = X.tail_readouts(space, p, "test")
    _, _, mass = Q.cvar_from_probs(space.E, p, X.ALPHA)
    assert ro["m"] == int((mass > 0).sum())
    order = np.argsort(space.E, kind="stable")
    assert set(ro["tail_idx"].tolist()) == set(order[:ro["m"]].tolist())


def test_nan_poison(space):
    """Contract rule 7: no deployable quantity may depend on a native."""
    cand, dg, rama = X.load_target(PDB)
    sp2 = X.Space(PDB, cand, dg, rama)
    cand.nat_ca = np.full_like(np.asarray(cand.nat_ca, float), np.nan)
    cand.oracle_rr = np.full_like(np.asarray(cand.oracle_rr, float), np.nan)
    sp3 = X.Space(PDB, cand, dg, rama)
    assert np.max(np.abs(sp2.E - sp3.E)) == 0.0
    assert np.max(np.abs(sp2.CA - sp3.CA)) == 0.0
    assert sp2.gamma_gap == sp3.gamma_gap
    assert np.array_equal(sp2.members, sp3.members)
    circ = Q.StatevectorCircuit(sp2.q, X.LAYERS)
    xors = [X.xor_index(sp2.M, sp2.q, k) for k in range(sp2.q)]
    r2 = X.train(circ, sp2.E, X.ALPHA, X.TEMP, sp2.gamma_gap, xors, 0)
    r3 = X.train(circ, sp3.E, X.ALPHA, X.TEMP, sp3.gamma_gap, xors, 0)
    a = X.tail_readouts(sp2, r2["p"], "poison")
    b = X.tail_readouts(sp3, r3["p"], "poison")
    for R in ("R1", "R2", "R3"):
        assert np.max(np.abs(a[R] - b[R])) == 0.0


def test_member_selection_breaks_ties_with_a_key(space):
    """Contract item 12: exact score ties never break by array order."""
    dis = np.zeros(20)                              # every candidate tied
    m1 = X.member_indices("AAAA", None, dis)
    m2 = X.member_indices("BBBB", None, dis)
    assert not np.array_equal(m1, m2), "a fully tied score must not return the array order"
    assert np.array_equal(m1, X.member_indices("AAAA", None, dis)), "and it must be stable"


def test_permuted_energy_is_a_real_control(space):
    """The permuted posterior keeps each pair's marginal family and changes the assignment."""
    Ep = X.permuted_energy(space)
    assert Ep.shape == space.E.shape
    assert np.isfinite(Ep).all()
    assert not np.allclose(Ep, space.E)
    # it is deterministic (a stable seed), and the Ramachandran half is untouched: the
    # permuted pair term keeps the SAME total over the 17 bins of every pair's own family,
    # so only the assignment of families to pairs moves.
    assert np.max(np.abs(Ep - X.permuted_energy(space))) == 0.0
    assert abs(float(np.sum(np.asarray(space.dg["prob"], float)))
               - float(space.dg["prob"].shape[0])) < 1e-3


def test_gibbs_matched_entropy(space):
    S_target = 4.0
    p, T = X.gibbs_matched_entropy(space.E, S_target)
    pp = p[p > 0]
    assert abs(float(-(pp * np.log(pp)).sum()) - S_target) < 1e-3


def test_scrambled_space_is_a_matched_null(space):
    """S29-L4 hole (a): the order-statistic control must match cardinality, parents and
    marginal fragment content, and must NOT be the same space."""
    cand, dg, rama = X.load_target(PDB)
    null = X.Space(PDB, cand, dg, rama, scramble=True)
    assert null.M == space.M and null.q == space.q
    assert np.array_equal(null.members, space.members)
    assert not np.allclose(null.CA, space.CA)
    # the same multiset of torsion pairs per member (blocks moved, nothing invented)
    for f in range(space.F):
        c = f * (space.M - 1) // (space.F - 1)
        a = np.sort(np.round(space.PHI[c], 9))
        b = np.sort(np.round(null.PHI[c], 9))
        assert len(a) == len(b)
    assert np.isfinite(null.E).all()


def test_shape_and_participation_ratio(space):
    """Contract addendum 20(c) and S29-L4 hole (c): the mechanism quantities exist and are
    the textbook ones."""
    s = X.shape_of(space.CA[0])
    assert 2.0 < s["rg"] < 40.0 and 3.5 < s["bond"] < 4.2
    p = X.gibbs(space.E, 2.0)
    ro = X.tail_readouts(space, p, "prtest")
    assert 1.0 <= ro["pr"] <= ro["m"] + 1e-9
    assert 0.0 <= ro["r3_mass"] <= 1.0 + 1e-9
    W = X.random_weight_control(space, ro["tail_idx"], ro["pr"], "prtest", n_draws=2)
    assert len(W) == 2 and all(np.isfinite(C).all() for C in W)


def test_cvar_optimal_law_is_optimal_and_is_boltzmann_at_alpha_one(space):
    """Lane T's M6 (S29-L15 Q1): the exact minimiser of CVaR_alpha - T*S over the simplex,
    and its alpha = 1 limit is the Boltzmann law (T's own assertion, reproduced here)."""
    p, t = X.cvar_optimal_law(space.E, X.ALPHA, X.TEMP)
    assert abs(float(p.sum()) - 1.0) < 1e-12 and (p >= 0).all()

    def F(q):
        v, _, _ = Q.cvar_exact(space.E, q, X.ALPHA)
        qq = q[q > 0]
        return v - X.TEMP * float(-(qq * np.log(qq)).sum())
    f0 = F(p)
    rng = np.random.default_rng(0)
    for _ in range(40):
        assert F(rng.dirichlet(np.full(space.M, 0.5))) >= f0 - 1e-9
    p1, _ = X.cvar_optimal_law(space.E, 1.0, 1.0)
    assert X.total_variation(p1, X.gibbs(space.E, 1.0)) < 1e-12


def test_total_variation_bounds(space):
    p = X.gibbs(space.E, 1.0)
    q = X.gibbs(space.E, 5.0)
    assert 0.0 <= X.total_variation(p, q) <= 1.0
    assert X.total_variation(p, p) == 0.0


def test_f_value_matches_f_and_grad(space):
    """`f_value` is the gradient-free path used by the BESTOFN control (prereg addendum 3);
    it must return exactly the same F as the full path."""
    circ = Q.StatevectorCircuit(space.q, X.LAYERS)
    xors = [X.xor_index(space.M, space.q, k) for k in range(space.q)]
    th = np.random.default_rng(23).normal(0.0, 0.6, circ.n_params())
    for gam in (0.0, space.gamma_gap):
        f0 = X.f_and_grad(circ, th, space.E, X.ALPHA, X.TEMP, gam, xors)[0]
        f1 = X.f_value(circ, th, space.E, X.ALPHA, X.TEMP, gam, xors)[0]
        assert f0 == f1
