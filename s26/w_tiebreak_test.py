#!/usr/bin/env python
"""s26/w_tiebreak_test.py -- synthetic tests for s26/w_tiebreak.py's draw logic.  Random
arrays only; no universe, no native."""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                   # noqa: E402

import w_tiebreak as T                               # noqa: E402


def synthetic_universe(nw=3000, k=500, seed=0):
    rng = np.random.default_rng(seed)
    sim = rng.integers(-10, 40, size=nw).astype(float)      # integer BLOSUM sums: many ties
    order = np.argsort(-sim, kind="stable")                  # the production convention
    return sim, order


def test_boundary_class_is_consistent():
    sim, order = synthetic_universe()
    prefix, tied, n_in = T.boundary_class(sim, order, k=500)
    s = sim[order[499]]
    assert (sim[prefix] > s).all() and len(prefix) + n_in == 500
    assert (sim[tied] == s).all() and n_in <= len(tied)
    # the production pool's tied members are the FIRST n_in of the tie class in corpus order
    assert np.array_equal(np.sort(order[:500][sim[order[:500]] == s]), tied[:n_in])
    print("  boundary class OK (tied %d, inside %d)" % (len(tied), n_in))


def test_draw_pool_is_a_valid_tiebreak_and_seeded():
    sim, order = synthetic_universe()
    pool, prefix, tied, n_in = T.draw_pool(sim, order, ("test", "X", 0), k=500)
    s = sim[order[499]]
    assert len(pool) == 500 and len(set(pool.tolist())) == 500
    assert set(prefix.tolist()) <= set(pool.tolist())
    drawn = pool[len(prefix):]
    assert (sim[drawn] == s).all() and len(drawn) == n_in
    pool2, *_ = T.draw_pool(sim, order, ("test", "X", 0), k=500)
    assert np.array_equal(pool, pool2)                   # same seed parts, same draw
    pool3, *_ = T.draw_pool(sim, order, ("test", "X", 1), k=500)
    assert not np.array_equal(pool, pool3)               # a different draw
    # every member of the tie class can appear: over 200 draws each tied member is drawn at least once
    seen = set()
    for kk in range(200):
        p, *_ = T.draw_pool(sim, order, ("test", "Y", kk), k=500)
        seen |= set(p[len(prefix):].tolist())
    assert seen == set(tied.tolist()), (len(seen), len(tied))
    print("  draw_pool OK")


if __name__ == "__main__":
    test_boundary_class_is_consistent()
    test_draw_pool_is_a_valid_tiebreak_and_seeded()
    print("ALL OK")
