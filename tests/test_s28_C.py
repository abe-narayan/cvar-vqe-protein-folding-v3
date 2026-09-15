"""tests/test_s28_C.py -- S28 lane C: the FAIL18 detector's feature functions and the
ranking-consuming readouts, on synthetic pools (no instrument, no native, no RMSD).

Covers: (1) features are finite and native-free (NaN-poison: `nat_ca` and `oracle_rr` set to
NaN give bit-identical features); (2) the readouts reduce to the uniform average at their
identity parameters and their weights are convex; (3) NaN-poison for every deployable readout
arm; (4) the nested logistic classifier finds a planted signal and not a random one; (5) the
permuted-ranker control is a different operator from the real one.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                 # noqa: E402
from s24 import d_harness as H                  # noqa: E402
from s27 import s28_C_fail18 as F               # noqa: E402
from s27 import s28_C_readout as R              # noqa: E402


def _synthetic_pool(seed=0, k=500, n=12):
    rng = np.random.default_rng(seed)
    base = np.cumsum(rng.normal(size=(n, 3)) * 2.5, 0)
    W = base[None] + rng.normal(scale=1.0, size=(k, n, 3))
    seq = "".join(rng.choice(list("ARNDCQEGHILKMFPSTWYV"), n))
    cand = H.Candidates(pdb="TEST", n=n, seq=seq, fold=0, W=W, nat_ca=base + 0.3,
                        oracle_rr=rng.random(k) * 3, source="synthetic")
    ch = {nm: rng.normal(size=k) for nm in ("DIS", "DISTPOT", "ENV", "CONTACT", "CONS")}
    ch["CONS"] = np.abs(ch["CONS"]) + 1.0
    sim = np.sort(rng.normal(size=k))[::-1] * 5 + 20
    npairs = n * (n - 1) // 2 - (n - 1)
    prob = rng.dirichlet(np.ones(17), size=npairs)
    dg = {"prob": prob}
    return cand, ch, sim, dg


def _poison(cand):
    return H.Candidates(pdb=cand.pdb, n=cand.n, seq=cand.seq, fold=cand.fold, W=cand.W.copy(),
                        nat_ca=np.full_like(cand.nat_ca, np.nan),
                        oracle_rr=np.full(cand.k, np.nan), source="poisoned")


# ------------------------------------------------------------------ Part 1: features
def test_features_finite_and_named():
    cand, ch, sim, dg = _synthetic_pool()
    f = F.features_one(cand, ch, sim, dg)
    assert set(f) == set(F.SP_NAMES) | set(F.CTRL_NAMES)
    assert all(np.isfinite(v) for v in f.values())
    assert 0.0 <= f["DISTPOT_overlap"] <= 1.0 and 0.0 <= f["sim_ent"] <= 1.0
    assert f["n"] == cand.n


def test_features_nan_poison_bit_identical():
    cand, ch, sim, dg = _synthetic_pool(seed=3)
    f0 = F.features_one(cand, ch, sim, dg)
    f1 = F.features_one(_poison(cand), ch, sim, dg)
    for k in f0:
        assert f0[k] == f1[k], k


def test_features_overlap_is_one_when_channel_equals_dis():
    cand, ch, sim, dg = _synthetic_pool(seed=4)
    ch = dict(ch)
    ch["DISTPOT"] = ch["DIS"].copy()
    f = F.features_one(cand, ch, sim, dg)
    assert f["DISTPOT_overlap"] == 1.0
    assert abs(f["DISTPOT_rho_dis"] - 1.0) < 1e-12


def test_auroc_and_logistic_on_planted_signal():
    rng = np.random.default_rng(0)
    n, p = 126, 6
    folds = np.repeat(np.arange(5), 26)[:n]
    X = rng.normal(size=(n, p))
    lin = X @ rng.normal(size=p)
    y = (lin > np.quantile(lin, 1 - 18 / 126)).astype(float)
    dec, thr, chosen = F.nested_logistic(X, y, folds, alphas=np.logspace(-1, 3, 5))
    assert F.auroc(y, dec) > 0.9
    assert len(chosen) == 5
    y_r = np.zeros(n)
    y_r[rng.choice(n, 18, replace=False)] = 1.0
    dec_r, _, _ = F.nested_logistic(X, y_r, folds, alphas=np.logspace(-1, 3, 5))
    assert F.auroc(y_r, dec_r) < 0.75
    assert F.auroc(np.array([1, 0, 1, 0]), np.array([0.9, 0.1, 0.8, 0.2])) == 1.0
    assert F.auroc(np.array([1, 0, 1, 0]), np.array([0.1, 0.9, 0.2, 0.8])) == 0.0


def test_single_rule_sign_is_nested():
    rng = np.random.default_rng(5)
    n = 126
    folds = np.repeat(np.arange(5), 26)[:n]
    x = rng.normal(size=n)
    y = (x < np.quantile(x, 18 / 126)).astype(float)      # LOW x = positive
    dec, thr = F.nested_single(x, y, folds)
    assert F.auroc(y, dec) > 0.95                          # the sign was learned on training folds


# ------------------------------------------------------------------ Part 2: readouts
def _set75(seed=1):
    rng = np.random.default_rng(seed)
    m, n = 75, 12
    base = np.cumsum(rng.normal(size=(n, 3)) * 2.0, 0)
    W = base[None] + rng.normal(scale=0.8, size=(m, n, 3))
    P = I.pairwise_rmsd(W)
    r = rng.random(m)
    return W, P, r


def test_identities_reduce_to_uniform_average():
    W, P, r = _set75()
    C0, _ = I.coordinate_average(W, P)
    Ca, idx = R.readout_mednb(W, P, r, 74)
    assert len(idx) == 75 and np.allclose(Ca, C0, atol=1e-9)
    Cb, keep = R.readout_trim(W, P, r, 0.0)
    assert len(keep) == 75 and np.allclose(Cb, C0, atol=1e-9)
    Cc, w = R.readout_divw(W, P, r, 0.0, 0.0)
    assert np.allclose(w, 1.0 / 75) and np.allclose(Cc, C0, atol=1e-9)


def test_readout_sizes_and_convexity():
    W, P, r = _set75(2)
    _, idx = R.readout_mednb(W, P, r, 20)
    assert len(idx) == 21 and idx[0] == int(np.argmin(r)) and len(set(idx.tolist())) == 21
    _, keep = R.readout_trim(W, P, r, 0.1)
    assert len(keep) == 67
    dropped = set(range(75)) - set(keep.tolist())
    assert dropped == set(np.argsort(-r)[:8].tolist())      # the 8 worst by the ranker
    for beta, gamma in ((1, 0), (0, 1), (2, 1)):
        w = R.divw_weights(P, r, beta, gamma)
        assert abs(w.sum() - 1) < 1e-12 and (w > 0).all()
    w = R.divw_weights(P, r, 1.0, 0.0)
    assert w[np.argmin(r)] > w[np.argmax(r)]                 # better ranker value, larger weight


def test_ranker_within_breaks_ties_by_key_not_order():
    cand, ch, sim, dg = _synthetic_pool(seed=7)
    top = np.arange(75)
    ch = dict(ch)
    ch["CONS"] = np.ones(cand.k)                              # all tied
    r1 = R.ranker_within(ch, top, "CONS", "A", perm=False)
    r2 = R.ranker_within(ch, top, "CONS", "B", perm=False)
    assert sorted(r1.tolist()) == list(range(75))
    assert not np.array_equal(r1, r2)                         # the key, not array order, decides
    assert not np.array_equal(r1, np.arange(75, dtype=float))


def test_readout_arms_nan_poison_bit_identical():
    cand, ch, sim, dg = _synthetic_pool(seed=11)
    pz = _poison(cand)
    for name, spec in R.all_arms():
        C0, i0 = R.emit(cand, ch, spec)
        C1, i1 = R.emit(pz, ch, spec)
        assert np.array_equal(C0, C1), name
        assert i0 == i1, name


def test_permuted_control_differs_from_real_operator():
    cand, ch, sim, dg = _synthetic_pool(seed=12)
    spec_real = ("trim", "CONS", False, 0.1)
    spec_perm = ("trim", "CONS", True, 0.1)
    C0, _ = R.emit(cand, ch, spec_real)
    C1, _ = R.emit(cand, ch, spec_perm)
    assert not np.allclose(C0, C1)
    Cp, _ = R.emit(cand, ch, None)
    Ci, _ = R.emit(cand, ch, ("divw", "CONS", False, (0.0, 0.0)))
    assert np.allclose(Cp, Ci, atol=1e-9)


def test_primary_chain_arms_exist():
    names = {a[0] for a in R.all_arms()}
    assert set(R.PRIMARY_CHAIN) <= names
    assert "PROD" in names and "DIVW[CONS,b=0,g=0]" in names
