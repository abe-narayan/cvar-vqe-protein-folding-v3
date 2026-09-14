#!/usr/bin/env python
"""s26/w_ensemble_test.py -- synthetic tests for s26/w_ensemble.py.  Random arrays only."""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                   # noqa: E402

import w_ensemble as E                               # noqa: E402
from s12 import instrument as I                      # noqa: E402
from core import data as D                           # noqa: E402


def test_blosum_matrices():
    M62 = E.blosum(62)
    assert np.array_equal(M62, np.asarray(D.BLOSUM62, float))
    for key in (45, 80):
        M = E.blosum(key)
        assert M.shape == (20, 20) and np.array_equal(M, M.T) and (np.diag(M) > 0).all()
    rng = np.random.default_rng(0)
    S = rng.integers(0, 20, size=(300, 11)).astype(np.int8); q = rng.integers(0, 20, size=11)
    assert np.array_equal(E.sims(62, q, S), D.blosum_similarity(q, S))
    print("  BLOSUM matrices OK")


def test_combine_is_frame_invariant_and_identity_on_copies():
    rng = np.random.default_rng(1)
    C = np.cumsum(rng.normal(size=(12, 3)) * 2.2, axis=0)
    # a rotated, translated copy combines to the reference itself
    th = 0.7; R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    C2 = C @ R.T + 5.0
    out = E.combine(C, [C2, C2])
    assert I.ca_rmsd(out, C) < 1e-9 and np.abs(out - C).max() < 1e-9
    # combining with a perturbed copy moves toward it by a third
    C3 = C + rng.normal(size=C.shape) * 0.5
    out = E.combine(C, [C3])
    assert I.ca_rmsd(out, C) < I.ca_rmsd(C3, C)
    print("  combine OK")


def test_union_and_boot_clouds():
    rng = np.random.default_rng(2)
    top = np.cumsum(rng.normal(size=(75, 10, 3)) * 2.2, axis=1)
    P = I.pairwise_rmsd(top); ref = top[I.medoid(P)]
    C = I.superpose_batch(top, ref).mean(0)
    U = E.union_cloud(np.concatenate([top, top, top], 0), ref)
    assert np.abs(U - C).max() < 1e-9                # a triple copy of the same set is the same cloud
    b1 = E.boot_cloud(top, np.random.default_rng(5)); b2 = E.boot_cloud(top, np.random.default_rng(5))
    assert np.allclose(b1, b2)                       # seeded
    b3 = E.boot_cloud(top, np.random.default_rng(6))
    assert not np.allclose(b1, b3)
    # a resample's cloud is nearer to the cloud than a typical member is (random-walk members here
    # are far more diverse than a real top-75, so the test is scale-relative, not an absolute A)
    member_dist = np.mean([I.ca_rmsd(w, C) for w in top])
    for b in (b1, b3):
        assert I.ca_rmsd(b, C) < member_dist, (I.ca_rmsd(b, C), member_dist)
    print("  union / boot OK (boot-cloud distance %.2f / %.2f against member distance %.2f)"
          % (I.ca_rmsd(b1, C), I.ca_rmsd(b3, C), member_dist))


if __name__ == "__main__":
    test_blosum_matrices()
    test_combine_is_frame_invariant_and_identity_on_copies()
    test_union_and_boot_clouds()
    print("ALL OK")
