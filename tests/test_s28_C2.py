"""tests/test_s28_C2.py -- S28 lane C2: the recognition audit's adapters and controls.

(1) On pool members the one-Context adapters reproduce the S27 pool-row channel values from
`s27/cache/<pdb>.npz` (skipped if the cache or the universe is absent); (2) the controls land
at the requested CA-RMSD; (3) the pairwise logistic finds a planted preference and not a
random one; (4) the Wilson interval and the fold-clustered CI behave.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s27 import s28_C2_recog_audit as R          # noqa: E402

PDB = "1A13"
CACHE = os.path.join(ROOT, "s27", "cache", PDB + ".npz")
UNIV = os.path.join(ROOT, "s8", "generate_univ", PDB + ".npz")
AMBER = os.path.join(ROOT, "s24", "cache_amber", PDB + ".npz")


@pytest.mark.skipif(not (os.path.exists(CACHE) and os.path.exists(UNIV) and os.path.exists(AMBER)),
                    reason="S27 channel cache / universe / AMBER cache not on disk")
def test_adapters_reproduce_pool_rows_on_pool_members():
    from s12 import instrument as I
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(PDB)
    u = I.load_univ(PDB)
    dg = I.distogram(PDB, cand.seq, cand.fold)
    idx = np.array([0, 7, 123])
    sc = R.ca_scores(cand, cand.W[idx], u, dg)
    for nm in ("DIS", "DIS_MEAN", "CONTACT_LL", "DISTPOT", "CONTACT", "ENV", "HP", "RG_LAW", "RG_UNIV",
               "EXVOL", "CAGEO", "SS_MATCH"):
        ref = np.asarray(ch[nm], float)[idx]
        assert np.allclose(sc[nm], ref, rtol=1e-9, atol=1e-9), nm
    k = cand.k
    assert np.allclose(sc["CONS_POOL"], np.asarray(ch["CONS"], float)[idx] * (k - 1) / k, rtol=1e-9, atol=1e-9)
    assert np.allclose(sc["DMAP_CONS_POOL"], np.asarray(ch["DMAP_CONS"], float)[idx], rtol=1e-9, atol=1e-9)
    assert np.allclose(sc["POOLGO_POOL"], np.asarray(ch["POOLGO"], float)[idx], rtol=1e-9, atol=1e-9)


@pytest.mark.skipif(not (os.path.exists(CACHE) and os.path.exists(UNIV) and os.path.exists(AMBER)),
                    reason="S27 channel cache / universe / AMBER cache not on disk")
def test_chain_adapters_reproduce_pool_rows_on_pool_members():
    from s12 import instrument as I
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(PDB)
    u = I.load_univ(PDB)
    dg = I.distogram(PDB, cand.seq, cand.fold)
    idx = np.array([3, 40])
    sc = R.chain_scores(cand, cand.W[idx], cand.PHI[idx], cand.PSI[idx], u, dg)
    for nm in ("RAMA", "DSSPHB", "ELEC", "TORS_CONS_POOL") + tuple("LEG_" + t for t in R.HL.LEG_TERMS):
        ref = np.asarray(ch[nm.replace("_POOL", "")], float)[idx]
        assert np.allclose(sc[nm], ref, rtol=1e-9, atol=1e-9), nm
    assert np.allclose(sc["LEG"], np.asarray(ch["LEG"], float)[idx], rtol=1e-6, atol=1e-6)


def _pair():
    rng = np.random.default_rng(3)
    n = 12
    base = np.cumsum(rng.normal(size=(n, 3)) * 2.5, 0)
    prod = base + rng.normal(scale=0.5, size=(n, 3))
    tgt = base + rng.normal(scale=3.0, size=(n, 3))
    return prod, tgt


def test_controls_land_at_the_requested_distance():
    prod, tgt = _pair()
    for d in (0.3, 1.0, 2.9):
        C, s = R.scale_to_distance(prod, tgt, d)
        assert abs(R.rmsd_between(C, prod) - d) < 1e-6 and s > 0
        Cg, _ = R.gauss_perturb(prod, d, "TEST", 1, "m")
        assert abs(R.rmsd_between(Cg, prod) - d) < 1e-6
    C0, s0 = R.scale_to_distance(prod, tgt, 0.0)
    assert np.allclose(C0, prod) and s0 < 1e-6


def test_gauss_draws_differ_and_are_stable():
    prod, _ = _pair()
    a, _ = R.gauss_perturb(prod, 0.3, "TEST", 0, "03")
    b, _ = R.gauss_perturb(prod, 0.3, "TEST", 1, "03")
    a2, _ = R.gauss_perturb(prod, 0.3, "TEST", 0, "03")
    assert not np.allclose(a, b) and np.array_equal(a, a2)


def test_pairwise_logistic_planted_vs_random():
    rng = np.random.default_rng(0)
    folds = np.repeat(np.arange(5), 26)[:126]
    D = rng.normal(size=(126, 5)); D[:, 2] -= 1.5
    dec, chosen = R.nested_pairwise(D, folds)
    assert float((dec < 0).mean()) > 0.85 and len(chosen) == 5
    Dr = rng.normal(size=(126, 5))
    dr, _ = R.nested_pairwise(Dr, folds)
    assert float((dr < 0).mean()) < 0.7
    # a GLOBAL relabelling (every pair flipped) is invisible to the sign-augmented model: the
    # learned w flips with it and the decision values are identical (the null must therefore
    # flip signs PER TARGET, which is what `analyse` does)
    dneg, _ = R.nested_pairwise(-D, folds)
    assert np.allclose(dneg, dec)
    eps = np.where(rng.random(126) < 0.5, 1.0, -1.0)
    dflip, _ = R.nested_pairwise(D * eps[:, None], folds)
    assert float((dflip < 0).mean()) < 0.7


def test_wilson_and_fold_ci():
    lo, hi = R.wilson(63, 126)
    assert lo < 0.5 < hi and abs((lo + hi) / 2 - 0.5) < 1e-9
    lo1, hi1 = R.wilson(126, 126)
    assert hi1 <= 1.0 and lo1 > 0.95
    folds = np.repeat(np.arange(5), 26)[:126]
    ind = np.ones(126)
    ci = R.fold_ci(ind, folds, n_boot=200)
    assert ci == [1.0, 1.0]
