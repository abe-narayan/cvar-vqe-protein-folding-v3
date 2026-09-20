"""tests/test_s29_O.py -- S29 lane O, the ORACLE ceiling ladder (`s29/s29_O_ladder.py`).

Synthetic pools only (no target, no native file, no cache): every rung reduces to the uniform
average at its identity parameter; the sparse NNLS reproduces the convex-hull rung at s = k;
the simplex solver matches a brute-force grid; the rigid-body removal annihilates rigid
displacements and leaves shape displacements alone; the leave-fold-out step never reads its
own fold and its emitted structure is bit-identical under the NaN-poison; the random-18 null.
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

from s12 import instrument as I                # noqa: E402
from s29 import s29_O_ladder as O              # noqa: E402


def _rot(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def _pool(k=40, n=9, seed=0, noise=0.8):
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    base = np.column_stack([2.3 * np.cos(np.radians(100) * t), 2.3 * np.sin(np.radians(100) * t), 1.5 * t])
    W = base[None] + rng.normal(0, noise, (k, n, 3))
    for a in range(k):
        W[a] = W[a] @ _rot(rng).T + rng.normal(0, 5, 3)
    nat = base + rng.normal(0, 0.5, (n, 3))
    return W, nat


# ------------------------------------------------------------------ identity reductions
def test_avg_of_equals_coordinate_average():
    W, nat = _pool()
    P = I.pairwise_rmsd(W)
    idx = np.array([3, 7, 11, 19, 20, 33])
    X, b = O.avg_of(W, idx, P)
    Xr, br = I.coordinate_average(W[idx])
    assert np.abs(X - Xr).max() < 1e-12 and idx[br] == b


def test_rung2_at_m_equals_k_is_the_uniform_average():
    W, nat = _pool(k=30)
    P = I.pairwise_rmsd(W)
    order = np.random.default_rng(1).permutation(len(W))
    curve = O.oracle_m_curve(W, order, P, nat, m_max=len(W))
    Xr, _ = I.coordinate_average(W)
    assert abs(curve[-1] - I.ca_rmsd(Xr, nat)) < 1e-12
    assert abs(curve[0] - I.ca_rmsd(W[order[0]], nat)) < 1e-12          # m = 1 is the member itself


def test_rung3_one_basin_is_the_uniform_average_and_k_basins_are_found():
    W, nat = _pool(k=30)
    P = I.pairwise_rmsd(W)
    r, X, size, nf, ne = O.oracle_best_basin(W, np.arange(len(W)), P, nat, 1)
    Xr, _ = I.coordinate_average(W)
    assert size == 30 and nf == 1 and abs(r - I.ca_rmsd(Xr, nat)) < 1e-12
    # two SHAPE-distinct families (a translation is invisible to Kabsch, so the second family is
    # a genuinely different conformation: the chain with its z axis reflected)
    rng = np.random.default_rng(4)
    Wb = W[:10].copy(); Wb[:, :, 2] *= -1.0
    W2 = np.concatenate([W[:20], Wb], 0)
    P2 = I.pairwise_rmsd(W2)
    lab = O.cluster_labels(P2, 2)
    assert len(set(lab.tolist())) == 2
    assert sorted(np.bincount(lab)[1:].tolist()) == [10, 20]
    # a basin of size 1 is not eligible: one random blob against ten near-identical chains
    W3 = np.concatenate([W[:10], rng.normal(0, 12, (1, W.shape[1], 3))], 0)
    P3 = I.pairwise_rmsd(W3)
    r, X, size, nf, ne = O.oracle_best_basin(W3, np.arange(11), P3, nat, 2)
    assert nf == 2 and ne == 1 and size == 10


def test_sparse_at_full_support_reproduces_the_hull():
    W, nat = _pool(k=12, n=8)
    from s27 import s28_A_amp as A
    frame = A.Frame(W, np.arange(len(W)))
    natp0 = O.superpose_one(nat, frame.ref)
    r_h, X_h, w_h = O.oracle_hull(frame.Wf, nat, 8, natp0=natp0, rounds=3)
    sp = O.oracle_sparse_greedy(frame.Wf, nat, 8, s_list=(12,), rounds=3, natp0=natp0)
    r_s, X_s, sup, w_s = sp[12]
    assert sorted(sup) == list(range(12))
    assert abs(r_s - r_h) < 1e-6
    # the greedy path is nested and monotone non-increasing in the FIXED-frame residual
    sp2 = O.oracle_sparse_greedy(frame.Wf, nat, 8, s_list=(2, 3, 5), rounds=0, natp0=natp0)
    assert set(sp2[2][2]) <= set(sp2[3][2]) <= set(sp2[5][2])


def test_simplex_solver_matches_brute_force_grid():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(9, 3)); y = A @ np.array([0.2, 0.5, 0.3]) + rng.normal(0, 0.3, 9)
    w = O.convex_nnls(A, y)
    assert abs(w.sum() - 1) < 1e-12 and (w >= 0).all()
    best = None
    g = np.linspace(0, 1, 201)
    for a in g:
        for b in g:
            if a + b <= 1 + 1e-12:
                v = np.linalg.norm(A @ np.array([a, b, 1 - a - b]) - y)
                if best is None or v < best:
                    best = v
    assert np.linalg.norm(A @ w - y) <= best + 1e-6


def test_hull_hits_a_target_inside_the_hull_and_uniform_at_identity():
    W, _ = _pool(k=10, n=8, seed=5)
    from s27 import s28_A_amp as A
    frame = A.Frame(W, np.arange(len(W)))
    w0 = np.random.default_rng(7).dirichlet(np.ones(10))
    target = (w0 @ frame.Wf).reshape(8, 3)
    r, X, w = O.oracle_hull(frame.Wf, target, 8, natp0=O.superpose_one(target, frame.ref))
    assert r < 1e-6
    # the uniform average is inside the hull: the solver reaches it too
    Xu = frame.Wf.mean(0).reshape(8, 3)
    r, X, w = O.oracle_hull(frame.Wf, Xu, 8, natp0=O.superpose_one(Xu, frame.ref))
    assert r < 1e-6


# ------------------------------------------------------------------ rung 6 geometry
def test_rigid_removal_annihilates_rigid_and_keeps_shape():
    W, nat = _pool(k=5, n=10)
    C = W[0]
    n = len(C)
    trans = np.tile([1.0, -2.0, 0.5], n)
    assert np.linalg.norm(O.remove_rigid(trans, C)) < 1e-10
    w = np.array([0.3, -0.2, 0.7])
    rot = np.cross(w[None, :], C - C.mean(0)).ravel()
    assert np.linalg.norm(O.remove_rigid(rot, C)) < 1e-10
    Q = O.rigid_basis(C)
    assert np.abs(Q.T @ Q - np.eye(6)).max() < 1e-10
    rng = np.random.default_rng(0)
    d = rng.normal(size=3 * n)
    d_shape = d - Q @ (Q.T @ d)
    assert np.abs(O.remove_rigid(d_shape, C) - d_shape).max() < 1e-10
    assert abs(O.cosine(d_shape, d_shape) - 1) < 1e-12
    # a rigid displacement contributes nothing to the cosine: adding one leaves cos unchanged
    assert abs(O.cosine(O.remove_rigid(d_shape + trans + rot, C), d_shape) - 1) < 1e-9


def test_axis_probe_u_is_native_free_and_t0_is_production():
    W, nat = _pool(k=20, n=9)
    C, _ = I.coordinate_average(W[:10])
    B, _ = I.coordinate_average(W[10:])
    rec, u = O.axis_probe(C, B, nat, "SYN", "LIB75")
    # u is a pure function of (C, B): recomputed here with no native in scope at all
    u_ref = O.remove_rigid((C - O.superpose_one(B, C)).ravel(), C)
    assert np.array_equal(u, u_ref)
    # and it does not change when the native does
    rec2, u2 = O.axis_probe(C, B, nat + 3.0, "SYN", "LIB75")
    assert np.array_equal(u, u2)
    assert abs(rec["r_t0"] - I.ca_rmsd(C, nat)) < 1e-12           # t = 0 is production
    Xm1 = C + (-1.0) * u.reshape(9, 3)
    assert abs(rec["curve"][0] - I.ca_rmsd(Xm1, nat)) < 1e-12
    assert -1 <= rec["cos_uv"] <= 1 and len(rec["curve"]) == len(O.T_GRID)


def test_lfo_choice_never_reads_its_own_fold_and_poison_is_bit_identical():
    rng = np.random.default_rng(11)
    n_t, folds = 25, np.repeat(np.arange(5), 5)
    T = O.T_GRID
    # every fold's curves have a different minimum; fold 0's own minimum is at t = 1.5 while
    # the other folds agree on t = 0.3, so the choice for fold 0 must be 0.3
    curves = np.empty((n_t, len(T)))
    for i in range(n_t):
        centre = 1.5 if folds[i] == 0 else 0.3
        curves[i] = (T - centre) ** 2 + rng.normal(0, 1e-3, len(T))
    t_of, t_fold, ties = O.lfo_choices(curves, folds, T)
    assert abs(t_fold[0] - 0.3) < 1e-9 and all(abs(t_of[folds == 0] - 0.3) < 1e-9)
    # poison: fold 0's curves NaN -> fold 0's choice unchanged (it never reads them)
    cp = curves.copy(); cp[folds == 0] = np.nan
    t_of2, t_fold2, _ = O.lfo_choices(cp, folds, T, only_folds=[0])
    assert t_fold2[0] == t_fold[0]
    # and a non-finite TRAINING slice is refused rather than silently producing a choice
    with pytest.raises(ValueError):
        O.lfo_choices(cp, folds, T)
    # the emitted structure is a pure function of (C, u, t)
    W, nat = _pool(k=6, n=7)
    C = W[0]; u = rng.normal(size=21)
    X1 = O.lfo_structure(C, u, t_fold[0]); X2 = O.lfo_structure(C, u, t_fold2[0])
    assert np.array_equal(X1, X2)


def test_rand18_null_detects_an_injected_stratum_effect():
    rng = np.random.default_rng(2)
    d = rng.normal(0, 0.1, 126)
    fm = np.zeros(126, bool); fm[:18] = True
    r0 = O.rand18_null(d, fm, n_draw=2000)
    assert 0.0 <= r0["p_one_sided"] <= 1.0 and r0["k"] == 18
    d2 = d.copy(); d2[fm] += 0.5
    r1 = O.rand18_null(d2, fm, n_draw=2000)
    assert r1["p_one_sided"] < 0.01 and r1["fail18_mean"] > r1["null_p97_5"]


def test_chain_groups_cover_every_item_once_and_rung6_first():
    assert len(O.ALL_ITEMS) == 33 and O.ALL_ITEMS[:3] == ["prod", "lfo_LIB75", "lfo_BPRIME"]
    assert len(set(O.ALL_ITEMS)) == len(O.ALL_ITEMS)
    assert set(O.CHAIN_GROUPS["E"]) == {"best1_top128", "bestm128", "hull_top128"}
    assert set(O.CHAIN_GROUPS["F"]) == {"hull_oa_top75", "hull_oa_pool"}


def test_hull_oa_is_at_least_as_good_as_the_shared_frame_hull_on_a_synthetic_pool():
    """`oa` poses every member on the target individually, so its reachable set contains the
    uniform average of the posed members and it can only be helped by the extra freedom; it is
    NOT emittable (posing needs the native) and the test records why the two rows differ."""
    W, nat = _pool(k=12, n=8, seed=9)
    from s27 import s28_A_amp as A
    frame = A.Frame(W, np.arange(len(W)))
    r_cf, _, _ = O.oracle_hull(frame.Wf, nat, 8, natp0=O.superpose_one(nat, frame.ref))
    r_oa, X, w = O.oracle_hull_oa(W, np.arange(len(W)), nat, 8)
    assert abs(w.sum() - 1) < 1e-12 and (w >= -1e-12).all()
    assert r_oa <= r_cf + 1e-9


def test_chain_shard_partition_is_exact_and_disjoint():
    pdbs = ["p%02d" % k for k in range(126)]
    seen = []
    for k in range(5):
        seen.append([p for j, p in enumerate(pdbs) if j % 5 == k])
    flat = [p for s_ in seen for p in s_]
    assert sorted(flat) == sorted(pdbs) and len(set(flat)) == 126
    for a in range(5):
        for b in range(a + 1, 5):
            assert not (set(seen[a]) & set(seen[b]))
    assert O.chain_rows_path(None) != O.chain_rows_path(0)
    assert O.chain_rows_path(3).endswith("_shard3.jsonl")
