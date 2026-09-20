"""Tests for s30 lane R -- the nativeness ladder instrument.

The three properties the verdict rests on:
  t1/t2  the SIZE-MATCHED twin is scale-invariant BY CONSTRUCTION (lane L's relay), and it makes
         the two PURE functions of Rg exactly constant;
  t3     the pool-referenced consistency channels REPRODUCE `ham_lib`'s own definitions when the
         candidate set IS the pool -- so the only thing the reimplementation changes is the
         reference index, which is the whole point;
  t4     no channel can read a native: the scoring namespace has no `nat_ca` / `oracle_rr` at all.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

PDB = "1A13"


@pytest.fixture(scope="module")
def ctx():
    from s24 import d_harness as H
    from s12 import instrument as I
    from s30 import s30_R_ladder as R
    cand = H.Candidates.from_universe(PDB)
    return dict(cand=cand, univ=I.load_univ(PDB),
                dg=I.distogram(PDB, cand.seq, int(cand.fold)),
                rcnt=R.rama_cnt()[int(cand.fold)], R=R)


def test_t1_si_twins_are_scale_invariant(ctx):
    """Uniformly rescaling every candidate must leave every *_SI channel unchanged, and must
    move at least one raw channel.  This is the property partialling cannot give you."""
    from types import SimpleNamespace
    R, cand = ctx["R"], ctx["cand"]
    base = SimpleNamespace(seq=cand.seq, n=int(cand.n), fold=int(cand.fold), k=int(cand.k),
                           W=np.asarray(cand.W, float), PHI=np.asarray(cand.PHI, float),
                           PSI=np.asarray(cand.PSI, float))
    blown = SimpleNamespace(**{**base.__dict__, "W": base.W * 1.137})     # pure scale, no shape
    a, _ = R.score_si(base, ctx["univ"], ctx["dg"], ctx["rcnt"])
    b, _ = R.score_si(blown, ctx["univ"], ctx["dg"], ctx["rcnt"])
    assert set(a) == set(b) and a
    raw_a = R.score_all(base, ctx["univ"], ctx["dg"], ctx["rcnt"])
    raw_b = R.score_all(blown, ctx["univ"], ctx["dg"], ctx["rcnt"])
    for k in a:
        # measured against the RAW channel's own scale: RG_UNIV_SI is identically ~0 after the
        # rescale (every candidate sits at the universe median), so its own sd is not a yardstick
        raw = raw_a[k[:-3]]
        scale = abs(np.nanmean(raw)) + np.nanstd(raw) + 1e-12
        d = np.nanmax(np.abs(a[k] - b[k])) / scale
        assert d < 1e-6, f"{k} is NOT scale invariant: {d}"
    moved = [k for k in raw_a if np.nanmax(np.abs(raw_a[k] - raw_b[k])) > 1e-6]
    assert "DISTPOT" in moved and "CONTACT" in moved, \
        "the raw pair potentials should move under a pure rescale -- that is lane L's point"


def test_t2_pure_rg_functions_become_constant(ctx):
    from types import SimpleNamespace
    R, cand = ctx["R"], ctx["cand"]
    base = SimpleNamespace(seq=cand.seq, n=int(cand.n), fold=int(cand.fold), k=int(cand.k),
                           W=np.asarray(cand.W, float), PHI=np.asarray(cand.PHI, float),
                           PSI=np.asarray(cand.PSI, float))
    si, _ = R.score_si(base, ctx["univ"], ctx["dg"], ctx["rcnt"])
    raw = R.score_all(base, ctx["univ"], ctx["dg"], ctx["rcnt"])
    for nm in ("RG_LAW", "RG_UNIV"):
        assert np.std(raw[nm]) > 1e-6, f"{nm} should vary over a real pool"
        rel = np.std(si[nm + "_SI"]) / (abs(np.mean(si[nm + "_SI"])) + 1e-12)
        assert rel < 1e-3, f"{nm}_SI is not constant (rel sd {rel:.2e})"


def test_t3_pool_reference_reproduces_ham_lib_on_the_pool(ctx):
    """With K = k the reimplementation must equal `ham_lib`'s own channels to float tolerance."""
    from s27 import ham_lib as HL
    R, cand = ctx["R"], ctx["cand"]
    cx = HL.Context(cand, ctx["univ"], ctx["dg"], ctx["rcnt"])
    mine = R.consistency_pool(cx.W, cx.PHI, cx.PSI, cx.D, cx.sep, cx.k)
    for nm in R.POOL_REF:
        ref = np.asarray(HL.CHANNELS[nm](cx), float)
        assert np.allclose(mine[nm], ref, rtol=1e-9, atol=1e-8), \
            f"{nm} differs from ham_lib: max {np.max(np.abs(mine[nm] - ref)):.3e}"


def test_t4_the_scoring_namespace_carries_no_native(ctx):
    """The object every channel is handed has no `nat_ca` and no `oracle_rr`, so a channel cannot
    read a native even by accident.  A construction argument, stronger than a poison run."""
    from types import SimpleNamespace
    from s27 import ham_lib as HL
    R, cand = ctx["R"], ctx["cand"]
    cl = SimpleNamespace(seq=cand.seq, n=int(cand.n), fold=int(cand.fold), k=int(cand.k),
                         W=np.asarray(cand.W, float), PHI=np.asarray(cand.PHI, float),
                         PSI=np.asarray(cand.PSI, float))
    assert not hasattr(cl, "nat_ca") and not hasattr(cl, "oracle_rr")
    cx = HL.Context(cl, ctx["univ"], ctx["dg"], ctx["rcnt"])
    assert not hasattr(cx, "nat_ca") and not hasattr(cx, "oracle_rr")
    # and the universe is consumed through exactly four keys, none of them the native
    u4 = {k: ctx["univ"][k] for k in ("W", "S", "PHI", "PSI")}
    cx2 = HL.Context(cl, u4, ctx["dg"], ctx["rcnt"])
    for nm, fn in HL.CHANNELS.items():
        v = np.asarray(fn(cx2), float)
        assert v.shape == (cl.k,), nm


def test_t5_ladder_rungs_have_ideal_geometry_and_graded_rmsd(ctx):
    """Every rung is built by the same function, so bond geometry is identical by construction;
    and one resampled residue already spans a wide range of RMSD-to-native."""
    from core import geometry as geo
    from s12 import instrument as I
    R, cand = ctx["R"], ctx["cand"]
    phi0, psi0, seq = R.native_torsions(PDB)
    s = R.Sampler(int(cand.fold), seq)
    rng = np.random.default_rng(0)
    PHI, PSI, MS = R.make_ladder(phi0, psi0, s, int(cand.n), rng, (1,), 64)
    W = np.asarray(geo.build_backbone_batch(PHI, PSI)["CA"], float)
    rr = np.asarray(I.kabsch_rmsd_batch(W, np.asarray(cand.nat_ca, float)), float)
    assert (MS == 1).all()
    assert rr.max() - rr.min() > 2.0, "a one-residue ladder should span a wide RMSD range"
    step = np.linalg.norm(np.diff(W, axis=1), axis=-1)
    assert step.std() < 0.05, "CA-CA steps must be identical across rungs (ideal geometry)"
