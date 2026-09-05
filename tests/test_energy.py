"""Equivalence tests for `core.energy` against the shipped Legacy stack.

The Legacy model is a *fitted* eleven-term potential. Unlike AMBER, nothing here can be
checked against an external reference implementation -- the only definition of "correct"
is the shipped one, so every term is compared against `energy_terms` structure by
structure, on real BLOSUM-pool geometry, and the measured difference is printed.

Assertions are `== 0.0`, not tolerances. `core.energy` lifts the terms byte-for-byte and
`components_batch` was deliberately built so its reductions keep the scalar path's order
(see that function's docstring for the two places where the obvious batching silently did
not). A tolerance here would hide exactly the reassociation this design exists to avoid.
"""
import os
import sys
import time

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import protein_geometry as geo                                       # noqa: E402
import peptide_db as db                                              # noqa: E402
import energy_terms as et                                            # noqa: E402
import legacy_field as lf                                            # noqa: E402

import core.energy as E                                              # noqa: E402

PID = "1A13"
K = 100


def _pool(pid=PID, k=K):
    from s5.lib import B62, encode, windows_full
    from s7.amber_native import pool_for
    p = db.by_pdb(pid)
    fold = db.folds(5)[p.seq]
    _, S, PHI, PSI, _ = windows_full(pool_for(fold, p.seq), p.n)
    sim = B62[S, encode(p.seq)[None, :]].sum(1)
    idx = np.argsort(-sim)[:k]
    phi, psi = PHI[idx].astype(float), PSI[idx].astype(float)
    return p, geo.build_backbone_batch(phi, psi), phi, psi


@pytest.fixture(scope="module")
def pool():
    return _pool()


def _singles(coords, B):
    return [{k: v[b] for k, v in coords.items()} for b in range(B)]


# ------------------------------------------------------------ the model itself
def test_eleven_terms_and_fitted_weights_are_intact():
    """The Legacy model is its eleven terms and its fitted weighting, or it is not it."""
    assert E.TERM_NAMES == et.TERM_NAMES
    assert len(E.TERM_NAMES) == 11
    assert E.DEFAULT_WEIGHTS == et.DEFAULT_WEIGHTS
    assert E.FITTED_WEIGHTS == lf.FITTED_WEIGHTS
    assert E.TERMS == lf.TERMS
    assert E.SUBSETS == lf.SUBSETS
    assert set(E.FREE_WEIGHTS) == set(et.FREE_WEIGHTS)
    assert set(E.COUPLED_TERMS) == set(et.COUPLED_TERMS)
    assert set(E.FILTER_TERMS) == set(et.FILTER_TERMS)
    print(f"\n11 terms: {E.TERM_NAMES}")
    print(f"FITTED_WEIGHTS: {E.FITTED_WEIGHTS}")


def test_sequence_tables_are_identical():
    """MJ contact potential, burial scale and formal charges, unchanged."""
    for seq in ("INWKGIAAMAKKLL", "ACDEFGHIKLMNPQRSTVWY", "GGGGGG"):
        b1, q1, m1 = et.sequence_arrays(seq, True)
        b2, q2, m2 = E.sequence_arrays(seq, True)
        assert np.array_equal(b1, b2) and np.array_equal(q1, q2)
        assert np.array_equal(m1, m2)
    assert et.MJ_CORRECTED == E.MJ_CORRECTED
    assert et.MJ_RAW == E.MJ_RAW
    assert et.BURIAL == E.BURIAL and et.CHARGE == E.CHARGE
    assert et.COULOMB == E.COULOMB and et.DIELECTRIC == E.DIELECTRIC
    assert et.BURIAL_NORM == E.BURIAL_NORM
    print(f"\nMJ table: {len(E.MJ_CORRECTED)} pairs, identical")


def test_scalar_energy_components_is_unchanged(pool):
    """The consolidated scalar path, term by term, on 100 real structures."""
    p, coords, phi, psi = pool
    singles = _singles(coords, K)
    worst = {}
    for b, s in enumerate(singles):
        o = et.energy_components(p.seq, s, phi=phi[b], psi=psi[b])
        n = E.energy_components(p.seq, s, phi=phi[b], psi=psi[b])
        assert list(o) == list(n)
        for t in E.TERM_NAMES:
            worst[t] = max(worst.get(t, 0.0), abs(o[t] - n[t]))
    print(f"\nscalar energy_components, {K} real structures:")
    for t in E.TERM_NAMES:
        print(f"  {t:18s} max |d| {worst[t]:.3e}")
    assert max(worst.values()) == 0.0


# ------------------------------------------------------------ the batched path
def test_components_batch_matches_scalar_exactly(pool):
    """All eleven terms, per structure, batched vs scalar. Must be 0.0, not small."""
    p, coords, phi, psi = pool
    singles = _singles(coords, K)
    old = [et.energy_components(p.seq, s, phi=phi[b], psi=psi[b])
           for b, s in enumerate(singles)]
    new = E.components_batch(p.seq, coords, phi=phi, psi=psi)

    assert list(new) == E.TERM_NAMES
    print(f"\ncomponents_batch vs energy_components, {K} real structures "
          f"({PID}, n={p.n}):")
    print(f"  {'term':18s} {'max abs diff':>13s} {'max rel diff':>13s}  differing")
    worst = 0.0
    for t in E.TERM_NAMES:
        o = np.array([c[t] for c in old], dtype=float)
        n = np.asarray(new[t], dtype=float)
        assert n.shape == (K,)
        da = np.abs(o - n)
        dr = da / np.maximum(np.abs(o), 1e-300)
        worst = max(worst, da.max())
        print(f"  {t:18s} {da.max():13.3e} {dr.max():13.3e}  {(o != n).sum():3d}/{K}")
    print(f"  {'ALL ELEVEN':18s} {worst:13.3e}")
    assert worst == 0.0

    # and the weighted combination that is what actually ranks structures
    wo = np.array([et.total_from_components(c, E.FITTED_WEIGHTS) for c in old])
    wn = E.totals_batch(new)
    print(f"  weighted total     max abs diff {np.abs(wo - wn).max():.3e}  "
          f"differing {(wo != wn).sum()}/{K}")
    assert np.array_equal(wo, wn)


def test_batched_dssp_matrix_matches_scalar(pool):
    """The batched H-bond geometry, against `protein_geometry` directly."""
    p, coords, _, _ = pool
    Eb, okb = E._dssp_matrices_batch(coords["N"], coords["C"], coords["O"])
    worst_e = 0.0
    bad_ok = 0
    for b in range(K):
        s = {k: coords[k][b] for k in ("N", "C", "O")}
        e0, ok0 = geo.dssp_energy_matrix(s, min_sep=2, cutoff=-0.5)
        worst_e = max(worst_e, float(np.abs(e0 - Eb[b]).max()))
        bad_ok += int(not np.array_equal(ok0, okb[b]))
    print(f"\nbatched DSSP matrix, {K} structures: max |dE| {worst_e:.3e}, "
          f"admissibility masks differing {bad_ok}")
    assert worst_e == 0.0
    assert bad_ok == 0


def test_batch_handles_missing_backbone_and_no_torsions(pool):
    """CA-only input must behave like the scalar path, not crash or invent terms."""
    p, coords, _, _ = pool
    ca_only = {"CA": coords["CA"][:4]}
    new = E.components_batch(p.seq, ca_only)
    for b in range(4):
        o = et.energy_components(p.seq, {"CA": coords["CA"][b]})
        for t in E.TERM_NAMES:
            assert o[t] == new[t][b], t
    assert np.all(new["torsion"] == 0.0)
    assert np.all(new["hbond_local"] == 0.0)
    print("\nCA-only pool: matches scalar on all 11 terms, torsion and hbond are 0")


def test_batch_rejects_unbatched_input(pool):
    p, coords, _, _ = pool
    with pytest.raises(ValueError):
        E.components_batch(p.seq, {"CA": coords["CA"][0]})
    with pytest.raises(ValueError):
        E.components_batch(p.seq[:-1], coords)


# ------------------------------------------------------------ generation field
def test_batch_legacy_field_is_carried_over(pool):
    """`BatchLegacy` is a generation field with documented form differences.

    It is NOT the scalar model and this test does not pretend it is: `steric` uses the
    same layout so it must match exactly, while `hbond` uses a different tie-break in the
    greedy match and `aromatic` uses the CB fallback. The point is that the consolidated
    copy behaves like the shipped one.
    """
    import torsion_lib2 as tl2
    p, coords, _, _ = pool
    tab = tl2.library_for(p.seq, 8, p.seq)
    rep_old = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
    rep_new = E.PerResidueTorsion(p.seq, E.library_for(p.seq, 8, p.seq),
                                  chi_bits=False)
    bo = lf.BatchLegacy(p.seq, rep_old)
    bn = E.BatchLegacy(p.seq, rep_new)
    to = bo.terms_from_coords(coords)
    tn = bn.terms_from_coords(coords)
    d = np.abs(to - tn).max()
    print(f"\nBatchLegacy.terms_from_coords: shape {tn.shape}, max |d| {d:.3e}")
    assert to.shape == tn.shape == (K, 11)
    assert d == 0.0

    # steric is the one term BatchLegacy shares with the scalar model by construction
    col = E.TERMS.index("steric")
    scalar = np.array([et.steric_term({k: v[b] for k, v in coords.items()}, p.seq)
                       for b in range(K)])
    print(f"  steric, batched field vs scalar term: "
          f"max |d| {np.abs(scalar - tn[:, col]).max():.3e}")


# ------------------------------------------------------------ speed
def test_components_batch_is_faster_than_the_scalar_loop(pool):
    """Report the speedup. The exactness tests above are what make it meaningful."""
    p, coords, phi, psi = pool
    singles = _singles(coords, K)

    def scalar():
        return [et.energy_components(p.seq, s, phi=phi[b], psi=psi[b])
                for b, s in enumerate(singles)]

    def batch():
        return E.components_batch(p.seq, coords, phi=phi, psi=psi)

    scalar(), batch()
    ts, tb = [], []
    for _ in range(5):
        t0 = time.perf_counter(); scalar(); ts.append(time.perf_counter() - t0)
        t0 = time.perf_counter(); batch(); tb.append(time.perf_counter() - t0)
    s, b = min(ts), min(tb)
    print(f"\n{PID} n={p.n}, {K} structures (best of 5):")
    print(f"  scalar loop  {s * 1e3:7.2f} ms  ({s / K * 1e3:.4f} ms/structure)")
    print(f"  batched      {b * 1e3:7.2f} ms  ({b / K * 1e3:.4f} ms/structure)")
    print(f"  speedup      {s / b:.2f}x")
    assert b < s
