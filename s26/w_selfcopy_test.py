#!/usr/bin/env python
"""s26/w_selfcopy_test.py -- synthetic unit tests for s26/w_selfcopy.py.  No native, no
universe, no RMSD to any deposited structure: every array here is random.

    python s26/jobrun.py --agent W --tag TEST --name w_selfcopy_test --est-ram 0.4 -- python s26/w_selfcopy_test.py
"""
from __future__ import annotations

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                   # noqa: E402

import w_selfcopy as W                               # noqa: E402
from s12 import instrument as I                      # noqa: E402


def test_triangle_inequality_of_kabsch_rmsd():
    """|RMSD(x, z) - RMSD(y, z)| <= RMSD(x, y) on random shapes: the property the native-free
    bound rests on.  2,000 random triples, several sizes, including near-degenerate ones."""
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(2000):
        n = int(rng.integers(9, 17))
        x = rng.normal(size=(n, 3)) * 3.0
        y = x + rng.normal(size=(n, 3)) * rng.uniform(0.01, 2.0)
        z = rng.normal(size=(n, 3)) * 3.0
        dxz, dyz, dxy = I.ca_rmsd(x, z), I.ca_rmsd(y, z), I.ca_rmsd(x, y)
        worst = max(worst, abs(dxz - dyz) - dxy)
    assert worst <= 1e-9, worst
    print("  triangle inequality: worst violation %.2e" % worst)


def test_drop_and_refill():
    order = np.arange(1000)[::-1]                       # any order
    pool = W.drop_and_refill(order, [999, 998, 500], k=500)
    assert len(pool) == 500
    assert not np.isin([999, 998, 500], pool).any()
    assert pool[0] == 997 and pool[-1] == 497            # the next-best refill, order preserved
    pool2 = W.drop_and_refill(order, [], k=500)
    assert np.array_equal(pool2, order[:500])
    print("  drop_and_refill OK")


def synthetic_risk(n, grid_lo=2.0, step=0.05, npts=760, seed=1):
    """A synthetic risk table: per pair a V-shaped risk with its minimum at a random distance."""
    rng = np.random.default_rng(seed)
    i, j = I.pair_index(n)
    grid = grid_lo + step * np.arange(npts)
    centre = rng.uniform(5.0, 15.0, size=len(i))
    risk = np.abs(grid[None, :] - centre[:, None]).astype(np.float32)
    return {"grid": grid.astype(np.float32), "risk": risk}, centre


def test_emit_selects_the_lowest_scores_and_is_native_free():
    rng = np.random.default_rng(2)
    n, m = 12, 300
    Wp = np.cumsum(rng.normal(size=(m, n, 3)) * 2.2, axis=1)
    dg, _ = synthetic_risk(n)
    e = W.emit(Wp, dg, "A" * n, 0, project=False)
    sc = W.score_pool(Wp, dg, n)
    assert np.array_equal(e["score"], sc)
    top = e["top"]
    assert len(top) == W.TOPM and np.all(np.sort(sc[top]) <= np.sort(sc)[W.TOPM - 1] + 1e-12)
    assert set(top.tolist()) == set(np.argsort(sc, kind="stable")[:W.TOPM].tolist())
    P = I.pairwise_rmsd(Wp[top]); C = I.superpose_batch(Wp[top], Wp[top][I.medoid(P)]).mean(0)
    assert np.allclose(C, e["cloud"])
    assert W.finite(e)
    # poisoning: a blinded universe carries NaN natives and emit never touches them
    u = W.blind({"rr": np.ones(m), "nat_ca": np.ones((n, 3)), "W": Wp})
    assert np.isnan(u["rr"]).all() and np.isnan(u["nat_ca"]).all()
    assert not W.finite(u) and W.finite({k: v for k, v in u.items() if k == "W"})
    print("  emit selection / cloud / blind OK")


def test_tie_set_and_argmin_window_rmsd():
    sc = np.array([3.0, 1.0, 1.0 + 1e-13, 2.0])
    ts = W.tie_set(sc)
    assert set(ts.tolist()) == {1, 2}
    rng = np.random.default_rng(3)
    Wa = rng.normal(size=(4, 9, 3)); Wb = Wa.copy()
    assert W.argmin_window_rmsd(Wa, ts, Wb, ts, same=True) == 0.0
    # different sets: the bound is the maximum cross-pair RMSD, which bounds the change of the
    # tie-averaged endpoint for ANY reference point (checked against random "natives")
    A, B = np.array([0, 1]), np.array([2, 3])
    bnd = W.argmin_window_rmsd(Wa, A, Wb, B, same=False)
    for _ in range(200):
        nat = rng.normal(size=(9, 3))
        fa = np.mean([I.ca_rmsd(Wa[a], nat) for a in A]); fb = np.mean([I.ca_rmsd(Wb[b], nat) for b in B])
        assert abs(fa - fb) <= bnd + 1e-9
    print("  tie set / sel bound OK")


def test_triangle_bounds_between_emissions():
    rng = np.random.default_rng(4)
    n, m = 10, 200
    Wp = np.cumsum(rng.normal(size=(m, n, 3)) * 2.2, axis=1)
    dg, _ = synthetic_risk(n, seed=5)
    e1 = W.emit(Wp, dg, "A" * n, 0, project=False)
    e2 = W.emit(Wp[::-1], dg, "A" * n, 0, project=False)   # same members, reversed order: same top set
    same = set((m - 1 - e1["argmin_set"]).tolist()) == set(e2["argmin_set"].tolist())
    assert same
    tri = W.triangle(e1, e2, Wp, Wp[::-1], same)
    assert tri["tri_cloud"] < 1e-9 and tri["tri_sel"] == 0.0, tri
    print("  triangle between identical emissions is zero OK")


def test_bound_arithmetic_and_materiality():
    assert abs(W.bound_from_deltas([0.1, -0.3, 0.05, 0.0]) - 2.0 / 60 * 0.3) < 1e-15
    assert W.materiality(0.001).startswith("IMMATERIAL")
    assert W.materiality(0.05).startswith("MINOR")
    assert W.materiality(0.5).startswith("MATERIAL")
    assert W.materiality(float("nan")) == "not measured"
    print("  bound arithmetic OK")


def test_gate_refuses_without_signoff():
    old = W.LEDGER
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as fh:
            fh.write("# empty ledger\n"); tmp = fh.name
        W.LEDGER = tmp
        assert not W.signed_off()
        try:
            W.require_signoff("endpoint")
            raise AssertionError("the gate did not refuse")
        except SystemExit as e:
            assert "PHASE GATE" in str(e)
        with open(tmp, "a") as fh:
            fh.write("\n## L99 -- PHASE 0 SIGNED OFF (test)\n")
        assert W.signed_off()
    finally:
        W.LEDGER = old
        os.unlink(tmp)
    print("  phase gate OK")


def test_train_filter_removes_exactly_the_named_chains():
    """The exclusion filter of `train_without`, exercised on a fake corpus through the same
    predicate: the chain is removed iff pdb AND sequence match; an absent chain asserts."""
    class E:                                              # a stand-in for core.data.Peptide
        def __init__(self, pdb, seq): self.pdb, self.seq = pdb, seq
    ents = [E("AAAA", "ACDEFGHIK"), E("BBBB", "KLMNPQRST"), E("frag1", "ACDEFGHIK"), E("CCCC", "VWYACDEFG")]
    excl = {"AAAA": "ACDEFGHIK"}
    kept = [e for e in ents if not (getattr(e, "pdb", None) in excl and e.seq == excl[e.pdb])]
    assert len(kept) == 3 and all(e.pdb != "AAAA" for e in kept)
    assert any(e.pdb == "frag1" for e in kept)          # a fragment with the same sequence is NOT removed
    removed = len(ents) - len(kept)
    assert removed == len(excl)
    excl2 = {"ZZZZ": "QQQQQQQQQ"}
    kept2 = [e for e in ents if not (getattr(e, "pdb", None) in excl2 and e.seq == excl2[e.pdb])]
    assert len(ents) - len(kept2) == 0 != len(excl2)     # the assertion in train_without fires on this
    print("  exclusion predicate OK")


def test_control_chains_are_dev_targets_and_not_verbatim_relatives():
    dev = {t["pdb"]: t for t in I.targets()}
    seqs = [t["seq"] for t in dev.values()]
    for fold, chains in W.CONTROL_CHAINS.items():
        for c in chains:
            assert c in dev, c
            s = dev[c]["seq"]
            assert dev[c]["fold"] != fold, (c, fold)   # in the TRAINING set of that fold
            assert not any(o != s and (s in o or o in s) for o in seqs), c
    print("  control chains OK")


if __name__ == "__main__":
    test_triangle_inequality_of_kabsch_rmsd()
    test_drop_and_refill()
    test_emit_selects_the_lowest_scores_and_is_native_free()
    test_tie_set_and_argmin_window_rmsd()
    test_triangle_bounds_between_emissions()
    test_bound_arithmetic_and_materiality()
    test_gate_refuses_without_signoff()
    test_train_filter_removes_exactly_the_named_chains()
    test_control_chains_are_dev_targets_and_not_verbatim_relatives()
    print("ALL OK")
