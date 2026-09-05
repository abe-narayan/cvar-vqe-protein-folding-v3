"""Equivalence and convention pins for `core.geometry`.

Every test here answers one of two questions:

  (a) does the consolidated implementation reproduce, ON REAL DATA, the one it replaces?
  (b) is a convention that has silently broken a result before now impossible to get wrong?

Tolerances are stated and justified individually. Most are exactly 0.0.

Run: ``python -m pytest tests/test_geometry.py -q``
"""
import glob
import math
import os
import random
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import geometry as g                                          # noqa: E402
from core import data as d                                             # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UNIV = os.path.join(ROOT, "s8", "generate_univ")
C2 = os.path.join(ROOT, "s8", "consensus2_cache")

#: The instrument the consolidated path must reproduce.
EXPECT_POOL_BEST = 1.7108
EXPECT_SHIPPED = 3.4540


def _universes(limit=None):
    if not os.path.isdir(UNIV):
        pytest.skip("s8/generate_univ absent (regenerable cache); nothing to check against")
    files = sorted(glob.glob(os.path.join(UNIV, "*.npz")))
    return files[:limit] if limit else files


@pytest.fixture(scope="module")
def peptides():
    return d.load()


# =====================================================================================
# 1. THE RMSD IMPLEMENTATION
# =====================================================================================
def _explicit_rotation_rmsd(W, ref):
    """The other form that was in use (`s8/audit8.rmsd_to`): build R, rotate, measure."""
    W = np.asarray(W, float)
    ref = np.asarray(ref, float)
    if W.ndim == 2:
        W = W[None]
    P = W - W.mean(1, keepdims=True)
    Q = ref - ref.mean(0, keepdims=True)
    U, S, Vt = np.linalg.svd(np.einsum("bni,nj->bij", P, Q))
    dd = np.sign(np.linalg.det(np.einsum("bji,bkj->bik", Vt, U)))
    E = np.zeros((len(W), 3, 3))
    E[:, 0, 0] = E[:, 1, 1] = 1.0
    E[:, 2, 2] = dd
    R = np.einsum("bji,bjk,blk->bil", Vt, E, U)
    return np.sqrt(((np.einsum("bij,bnj->bni", R, P) - Q[None]) ** 2).sum(-1).mean(-1))


def test_ca_rmsd_batch_is_bit_identical_to_protein_geometry():
    """core.geometry.ca_rmsd_batch vs the module it replaces, on real pool windows.

    Not "close": exactly equal. Both are the singular-value residual with the same
    operation order, so any difference at all would mean one of them had been changed.
    """
    import protein_geometry as pg
    worst = 0.0
    n = 0
    for f in _universes(40):
        z = np.load(f, allow_pickle=True)
        W = np.asarray(z["W"], float)[np.asarray(z["order"])[:400]]
        nat = np.asarray(z["nat_ca"], float)
        worst = max(worst, float(np.abs(g.ca_rmsd_batch(W, nat)
                                        - pg.ca_rmsd_batch(W, nat)).max()))
        n += len(W)
    assert n > 10000
    assert worst == 0.0, f"differs on {n} superpositions by {worst}"


def test_ca_rmsd_batch_reproduces_the_cached_rr_of_every_universe():
    """The 126-target instrument's own numbers, recomputed from coordinates.

    ``rr`` is stored as float32, so the tolerance is float32 round-off on a ~2 A value
    (6e-07), not a modelling tolerance. Measured worst over all 126 targets: 4.8e-07 --
    which is the number an earlier audit reported as a *disagreement between Kabsch
    implementations*. It is not; it is the cache's dtype.
    """
    worst = 0.0
    for f in _universes():
        z = np.load(f, allow_pickle=True)
        o = np.asarray(z["order"])[:500]
        mine = g.ca_rmsd_batch(np.asarray(z["W"], float)[o], np.asarray(z["nat_ca"], float))
        worst = max(worst, float(np.abs(mine - np.asarray(z["rr"], float)[o]).max()))
    assert worst < 1e-6, worst


def test_singular_value_and_explicit_rotation_forms_agree():
    """The two Kabsch forms in the repository, on real data.

    8.9e-14 over 8.0 M superpositions -- double-precision round-off through an SVD, and
    the justification for keeping the cheaper singular-value form.
    """
    worst = 0.0
    for f in _universes(30):
        z = np.load(f, allow_pickle=True)
        W = np.asarray(z["W"], float)[np.asarray(z["order"])[:300]]
        nat = np.asarray(z["nat_ca"], float)
        worst = max(worst, float(np.abs(g.ca_rmsd_batch(W, nat)
                                        - _explicit_rotation_rmsd(W, nat)).max()))
    assert worst < 1e-11, worst


def test_pool_best_and_shipped_selected_reproduce_the_instrument():
    """1.7108 and 3.4540, computed through the consolidated geometry."""
    if not os.path.isdir(C2):
        pytest.skip("s8/consensus2_cache absent")
    best, sel = [], []
    for f in sorted(glob.glob(os.path.join(C2, "*.npz"))):
        z = np.load(f, allow_pickle=True)
        rr = g.ca_rmsd_batch(np.asarray(z["W"], float), np.asarray(z["nat_ca"], float))
        best.append(rr.min())
        sel.append(rr[int(np.argmin(np.asarray(z["sc"], float)))])
    assert abs(float(np.mean(best)) - EXPECT_POOL_BEST) < 1e-3, np.mean(best)
    assert abs(float(np.mean(sel)) - EXPECT_SHIPPED) < 1e-3, np.mean(sel)


def test_ca_rmsd_forbids_reflection():
    """A mirror image must score as different. This is the whole reason the sign flip on
    the smallest singular value is there: a distance matrix cannot tell them apart, and
    15-24% of some candidate pools are left-handed copies."""
    rng = np.random.default_rng(0)
    for p in d.load()[:20]:
        if p.n < 6:
            continue
        m = p.ca.copy()
        m[:, 2] *= -1.0
        assert float(g.ca_rmsd_batch(m[None], p.ca)[0]) > 0.5
    # and a pure rotation must score zero, to the residual form's noise floor
    q = d.load()[0].ca
    R = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    if np.linalg.det(R) < 0:
        R[:, 0] *= -1
    assert float(g.ca_rmsd_batch((q @ R.T + 7.0)[None], q)[0]) < 1e-6


def test_the_residual_forms_noise_floor_at_zero_is_known_and_bounded():
    """A property of the form that is kept, stated rather than discovered later.

    ``|P|^2 + |Q|^2 - 2 sum(s)`` is a difference of two quantities of order ``n * r^2``,
    so at RMSD 0 it cancels to about ``sqrt(eps) * r`` -- measured at 1.9e-07 A on a real
    13-mer, against 1e-15 for the explicit-rotation form, which never subtracts large
    numbers. That is nine orders below the 1-4 A scale every result in this repository is
    reported at, and it buys a third of the flops. It matters only if someone ever tries to
    use this to detect exact duplicates, so: DON'T -- compare coordinates for that.
    """
    rng = np.random.default_rng(11)
    q = d.load()[0].ca
    R = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    if np.linalg.det(R) < 0:
        R[:, 0] *= -1
    P = (q @ R.T + 7.0)[None]                       # the same structure, rigidly moved
    same = float(g.ca_rmsd_batch(P, q)[0])
    exact = float(_explicit_rotation_rmsd(P, q)[0])
    assert same < 1e-6, same                        # the floor, and it is bounded
    assert exact < 1e-12, exact                     # the form that never cancels
    # Away from zero the two forms agree to double-precision round-off, which is the
    # actual claim: the floor exists only where the answer is itself ~0.
    W = np.stack([p.ca for p in d.load()[:8] if p.n == len(q)][:4])
    if len(W) >= 2:
        assert np.abs(g.ca_rmsd_batch(W, q)
                      - _explicit_rotation_rmsd(W, q)).max() < 1e-11


def test_weighted_rmsd_reduces_to_the_unweighted_one():
    """At uniform weights the weighted form must BE the unweighted one -- the only thing
    that makes a re-weighting bound comparable with the number it is supposed to bound."""
    for p in d.load()[:15]:
        W = np.stack([p.ca[::-1], np.roll(p.ca, 2, axis=0)])
        w = np.ones(p.n)
        # Away from zero: double-precision agreement, the claim that matters.
        assert np.abs(g.wca_rmsd_batch(W, p.ca, w)
                      - g.ca_rmsd_batch(W, p.ca)).max() < 1e-11
        # At zero both forms sit on the cancellation floor characterised above.
        assert np.abs(g.wca_rmsd_batch(p.ca[None], p.ca, w)
                      - g.ca_rmsd_batch(p.ca[None], p.ca)).max() < 1e-6


def test_rotation_onto_matches_kabsch_superpose_side():
    """THE TRANSPOSED-ROTATION PIN.

    `kabsch_superpose` applies its rotation on the left; `rotation_onto` returns the matrix
    to POST-multiply by. Post-multiplying by the wrong one is still a valid rotation, so
    nothing raises -- it silently solves everything downstream in a mis-rotated frame,
    which is how a weighting solve returned a bound worse than the equal-weight mean it
    relaxes (2.648 -> 2.087 once fixed). If this test fails, that class of bug is back.
    """
    rng = np.random.default_rng(3)
    for p in d.load()[:12]:
        X = p.ca + rng.normal(0, 0.4, p.ca.shape)
        Y = p.ca
        Xc, Yc = X - X.mean(0), Y - Y.mean(0)
        M = g.rotation_onto(Xc, Yc)
        by_convention = Xc @ M + Y.mean(0)
        by_superpose = g.kabsch_superpose(X, Y)
        assert np.abs(by_convention - by_superpose).max() < 1e-9
        # and the transpose is measurably WORSE, i.e. the test has teeth
        wrong = Xc @ M.T + Y.mean(0)
        assert g.rmsd(wrong, Y) >= g.rmsd(by_superpose, Y) - 1e-12


def test_pairwise_ca_rmsd_is_symmetric_with_a_zero_diagonal():
    W = np.stack([p.ca for p in d.load()[:24] if p.n == d.load()[0].n][:8])
    if len(W) < 3:
        pytest.skip("not enough equal-length peptides")
    M = g.pairwise_ca_rmsd(W)
    assert np.abs(M - M.T).max() < 1e-12
    assert np.abs(np.diag(M)).max() < 1e-6      # residual-form floor, see above


# =====================================================================================
# 2. TORSIONS AND BACKBONE
# =====================================================================================
def test_dihedral_batch_is_exact_against_the_scalar_form(peptides):
    """Difference exactly 0.0 over every backbone torsion in the database."""
    worst = 0.0
    for p in peptides[:120]:
        b = g.build_backbone(p.phi, p.psi)
        N, CA, C = b["N"], b["CA"], b["C"]
        batch = g.dihedral_batch(C[:-1], N[1:], CA[1:], C[1:])
        scalar = np.array([g.dihedral(C[i - 1], N[i], CA[i], C[i])
                           for i in range(1, len(CA))])
        worst = max(worst, float(np.abs(batch - scalar).max()))
    assert worst == 0.0, worst


def test_extract_torsions_matches_the_legacy_loop(peptides):
    """The vectorised extractor against `protein_geometry`'s per-residue loop."""
    import protein_geometry as pg
    worst = 0.0
    for p in peptides[:200]:
        b = g.build_backbone(p.phi, p.psi)
        a1, a2 = pg.extract_torsions(b["N"], b["CA"], b["C"])
        b1, b2 = g.extract_torsions(b["N"], b["CA"], b["C"])
        worst = max(worst, float(np.abs(a1 - b1).max()), float(np.abs(a2 - b2).max()))
    assert worst == 0.0, worst


def test_extract_torsions_batch_matches_the_single_form(peptides):
    ps = [p for p in peptides if p.n == 13][:32]
    if len(ps) < 4:
        pytest.skip("not enough length-13 peptides")
    bb = [g.build_backbone(p.phi, p.psi) for p in ps]
    N = np.stack([b["N"] for b in bb])
    CA = np.stack([b["CA"] for b in bb])
    C = np.stack([b["C"] for b in bb])
    phi, psi = g.extract_torsions_batch(N, CA, C)
    for k in range(len(ps)):
        a, b = g.extract_torsions(N[k], CA[k], C[k])
        assert np.abs(phi[k] - a).max() == 0.0
        assert np.abs(psi[k] - b).max() == 0.0


def test_batched_backbone_matches_the_scalar_builder(peptides):
    """Not approximately -- exactly. The search evaluates structures in batches and the
    batched builder is the hot path; a drift here changes every energy."""
    ps = [p for p in peptides if p.n == 12][:48]
    if len(ps) < 4:
        pytest.skip("not enough length-12 peptides")
    phi = np.stack([p.phi for p in ps])
    psi = np.stack([p.psi for p in ps])
    B = g.build_backbone_batch(phi, psi)
    for k, p in enumerate(ps):
        s = g.build_backbone(p.phi, p.psi)
        for key in ("N", "CA", "C", "CB", "O"):
            assert np.abs(B[key][k] - s[key]).max() == 0.0, key


def test_torsion_roundtrip_is_stable(peptides):
    """Build from torsions, read the torsions back: the interior ones must return."""
    for p in peptides[:40]:
        b = g.build_backbone(p.phi, p.psi)
        phi, psi = g.extract_torsions(b["N"], b["CA"], b["C"])
        assert np.abs(np.angle(np.exp(1j * (phi[1:] - p.phi[1:])))).max() < 1e-8
        assert np.abs(np.angle(np.exp(1j * (psi[:-1] - p.psi[:-1])))).max() < 1e-8


# =====================================================================================
# 3. PDB IO
# =====================================================================================
def _sample_pdbs(k=120, seed=0):
    paths = (sorted(glob.glob(os.path.join(ROOT, "pdbs", "*.pdb")))
             + sorted(glob.glob(os.path.join(ROOT, "pdbs_ext", "*.pdb")))
             + sorted(glob.glob(os.path.join(ROOT, "prots", "*.pdb"))))
    if not paths:
        pytest.skip("no PDB corpora present")
    random.Random(seed).shuffle(paths)
    return paths[:k]


def test_fast_pdb_reader_matches_biopython_bit_for_bit():
    """The column reader against `Bio.PDB.PDBParser` on all three corpora.

    Exact on 61/61 `pdbs/`, 1462/1463 `pdbs_ext/` and 600/600 sampled `prots/`; the single
    outlier returns no standard residues under BOTH readers, so it is agreement, not a
    mismatch. The two disorder conventions this has to get right are Biopython's own:
    a disordered ATOM takes the first altloc of the highest occupancy, and a disordered
    RESIDUE (1CBN's SER/PRO at position 22) takes the LAST residue name deposited.
    """
    checked = 0
    for path in _sample_pdbs(120):
        try:
            ref = g._read_models_biopython(path)
        except Exception:
            continue
        mine = g._read_models(path)

        def trim(ms):
            out, r = [], None
            for m in ms:
                if r is None:
                    r = m[0]
                elif m[0] != r:
                    break
                out.append(m)
            return out

        ref, mine = trim(ref), trim(mine)
        assert len(ref) == len(mine), path
        for a, b in zip(ref, mine):
            assert a[0] == b[0], f"{path}: sequence differs"
            for k in (1, 2, 3):
                assert np.abs(a[k] - b[k]).max() == 0.0, f"{path}: coordinates differ"
        checked += 1
    assert checked > 50


def test_pdb_cache_returns_what_was_parsed(tmp_path, monkeypatch):
    """A cache that returned something else would be invisible, so it is checked."""
    from core import cache as C
    monkeypatch.setattr(C, "ROOT", str(tmp_path))
    C._MEM.clear()
    del C._MEM_ORDER[:]
    path = _sample_pdbs(1)[0]
    cold = g.parse_pdb_ensemble(path)
    C._MEM.clear()
    del C._MEM_ORDER[:]
    warm = g.parse_pdb_ensemble(path)
    assert len(cold) == len(warm)
    for a, b in zip(cold, warm):
        assert a[0] == b[0]
        for k in (1, 2, 3):
            assert np.abs(a[k] - b[k]).max() == 0.0


def test_pdb_access_is_logged():
    """Validation proves the optimiser never touches a native by reading this log."""
    g.reset_pdb_log()
    assert g.get_pdb_log() == []
    g.parse_pdb_ensemble(_sample_pdbs(1)[0])
    assert len(g.get_pdb_log()) == 1


# =====================================================================================
# 4. CONTACTS, SECONDARY STRUCTURE, CA-TRACE PRIOR, FLOOR
# =====================================================================================
def test_contacts_and_secondary_structure_match_the_legacy_module(peptides):
    import protein_geometry as pg
    for p in peptides[:60]:
        c = g.build_backbone(p.phi, p.psi)
        assert g.contact_map(c["CA"]) == pg.contact_map(c["CA"])
        assert g.assign_secondary_structure(c) == pg.assign_secondary_structure(c)
        E1, ok1 = g.dssp_energy_matrix(c)
        E2, ok2 = pg.dssp_energy_matrix(c)
        # Residue 0 has no amide H, so the first row is NaN by construction in both.
        assert np.array_equal(E1, E2, equal_nan=True)
        assert (ok1 == ok2).all()


def test_ss_agreement_matches_the_scalar_definition():
    a, b = "HHHECCCEEH", "HHCECCCEEE"
    assert abs(g.ss_agreement(a, b)
               - sum(1 for x, y in zip(a, b) if x == y) / len(a)) < 1e-12
    assert g.ss_agreement("", "HH") == 0.0


def test_pseudo_angles_match_catrace(peptides):
    import catrace
    for p in peptides[:60]:
        t1, u1 = g.pseudo_angles(p.ca)
        t2, u2 = catrace.pseudo_angles(p.ca[None])
        assert np.abs(t1 - t2).max() == 0.0
        assert np.abs(u1 - u2).max() == 0.0


def test_catrace_prior_scores_a_mirror_differently(peptides):
    """The term exists to break the distance matrix's reflection blind spot; if the mirror
    penalty were zero it would not do that, so it is measured rather than asserted."""
    prior = g.CATracePrior()
    ca = np.stack([p.ca for p in peptides if p.n == 13][:16])
    if len(ca) < 4:
        pytest.skip("not enough length-13 peptides")
    pen = prior.mirror_penalty(ca)
    assert np.abs(pen).mean() > 1e-3
    hand = prior.handedness(ca)
    mir = ca.copy()
    mir[..., 2] *= -1.0
    assert np.abs(prior.handedness(mir) + hand).max() < 1e-12


def test_representation_floor_matches_the_legacy_module(peptides):
    import floor as legacy_floor
    p = [q for q in peptides if 9 <= q.n <= 11][0]
    rep = d.make_representation("torsion", p.n, 4, sequence=p.seq)
    a = g.representation_floor(rep, p.ca, p.phi, p.psi, restarts=2, seed=0)
    b = legacy_floor.floor(rep, p.ca, p.phi, p.psi, restarts=2, seed=0)
    assert abs(a["projection"] - b["projection"]) < 1e-12
    assert abs(a["floor"] - b["floor"]) < 1e-12
    assert (a["states"] == b["states"]).all()


def test_project_states_matches_the_legacy_loop(peptides):
    import floor as legacy_floor
    for p in peptides[:40]:
        rep = d.make_representation("torsion", p.n, 4, sequence=p.seq)
        assert (g.project_states(rep, p.phi, p.psi)
                == legacy_floor.project(rep, p.phi, p.psi)).all()


# =====================================================================================
# 5. THE numerics CONTRACT
# =====================================================================================
def test_numerics_contract_is_the_canonical_implementation():
    """`core/__init__.py` routes the `numerics` backend here. The aliases must be the same
    objects `s7/audit.py` retyped by hand, or the consolidation is cosmetic."""
    import s7.audit as audit
    z = np.load(_universes(1)[0], allow_pickle=True)
    W = np.asarray(z["W"], float)[:200]
    nat = np.asarray(z["nat_ca"], float)
    assert np.abs(g.kabsch_rmsd_batch(W, nat) - audit.kabsch_rmsd_batch(W, nat)).max() == 0.0
    assert np.abs(g.B62 - audit.B62).max() == 0.0
    seq = str(z["seq"])
    assert (g.encode(seq) == audit.encode(seq)).all()
    i1, j1 = g.pair_index(len(seq), 2)
    i2, j2 = audit.pair_index(len(seq), 2)
    assert (i1 == i2).all() and (j1 == j2).all()
    assert np.abs(g.pair_dists(W, i1, j1) - audit.pair_dists(W, i2, j2)).max() == 0.0
