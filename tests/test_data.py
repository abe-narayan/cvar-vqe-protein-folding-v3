"""Equivalence and convention pins for `core.data` and `core.cache`.

The four conventions pinned here are the four that have each silently produced a wrong
result in this repository:

  1. the alphabet. Two 20-letter orders were in use and they agree on three letters. An
     audit decoded 126 targets' cached window banks with the wrong one and measured every
     identity against a scrambled sequence. Pinned against the deposited banks themselves.
  2. THE IDENTITY NORMALISATION. `identity` divides by the LONGER sequence, which is leaky
     at the member level -- 1CEK's 13 residues sit verbatim inside the 25-residue 1A11 and
     score 0.520. Results on record depend on this, so it is preserved and pinned, and the
     containment form is pinned separately so the two can never be swapped.
  3. RETRIEVAL TIE-BREAKING. Stable argsort reproduces the deposited pools; unstable moves
     up to 47 of 500 members.
  4. THE CACHE KEY. A key that omits a parameter is a scientifically invalid collision and
     is silent. `CacheCollision` makes it loud.

Run: ``python -m pytest tests/test_data.py -q``
"""
import glob
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import cache as C                                            # noqa: E402
from core import data as d                                             # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UNIV = os.path.join(ROOT, "s8", "generate_univ")


def _universes(limit=None):
    if not os.path.isdir(UNIV):
        pytest.skip("s8/generate_univ absent (regenerable cache)")
    files = sorted(glob.glob(os.path.join(UNIV, "*.npz")))
    return files[:limit] if limit else files


@pytest.fixture(scope="module")
def peptides():
    return d.load()


# 1. the alphabet
def test_the_canonical_alphabet_is_the_one_the_banks_are_written_in():
    """THE PIN. Decode cached window banks and require the result to be real sequence.

    Every window in a target's universe is a contiguous slice of some library peptide, so a
    correctly decoded window is a real peptide substring. Under `ALPHABET` the exact-sequence
    windows decode to their own target; under `ALPHABET_ALT` they decode to noise. This is
    the check the void audit did not have.
    """
    ok_canon = ok_alt = 0
    checked = 0
    for f in _universes():
        z = np.load(f, allow_pickle=True)
        seq = str(z["seq"])
        S = np.asarray(z["S"])
        top = np.asarray(z["order"])[:200]
        # the single most BLOSUM-similar window is, for a leaked target, the target itself
        if any(d.decode(S[k], d.ALPHABET) == seq for k in top):
            ok_canon += 1
        if any(d.decode(S[k], d.ALPHABET_ALT) == seq for k in top):
            ok_alt += 1
        checked += 1
    assert checked >= 100
    # Four of the 126 targets have their own sequence verbatim in the pool. Under the wrong
    # alphabet, none of them do -- which is exactly how the wrong audit read 0/126.
    assert ok_canon >= 4, ok_canon
    assert ok_alt == 0, ok_alt


def test_the_two_alphabets_agree_on_exactly_three_letters():
    """A, S and T. That is why the failure was silent: most letters still decoded to a
    letter, and only the identities were wrong."""
    agree = [a for a, b in zip(d.ALPHABET, d.ALPHABET_ALT) if a == b]
    assert agree == ["A", "S", "T"], agree
    assert sorted(d.ALPHABET) == sorted(d.ALPHABET_ALT)


def test_the_1cek_case_reproduces_exactly():
    """The documented instance, to three decimals, so this file records the evidence.

    1CEK's own exact-sequence pool window decodes to AISVLLAQAVFLL in the bank's alphabet
    (positional identity 1.000) and to ALSYMMAGAYQMM in the other (0.308).
    """
    f = os.path.join(UNIV, "1CEK.npz")
    if not os.path.exists(f):
        pytest.skip("1CEK universe absent")
    z = np.load(f, allow_pickle=True)
    seq = str(z["seq"])
    S = np.asarray(z["S"])
    hits = [k for k in range(len(S)) if d.decode(S[k], d.ALPHABET) == seq]
    assert len(hits) == 1
    k = hits[0]
    assert d.decode(S[k], d.ALPHABET) == "AISVLLAQAVFLL"
    assert d.decode(S[k], d.ALPHABET_ALT) == "ALSYMMAGAYQMM"
    frac = sum(a == b for a, b in zip(seq, d.decode(S[k], d.ALPHABET_ALT))) / len(seq)
    assert abs(frac - 0.3077) < 5e-4, frac


def test_decode_cannot_be_called_without_naming_an_alphabet():
    """The structural defence: the failure needed a default, so there is no default."""
    with pytest.raises(TypeError):
        d.decode(np.array([0, 1, 2]))          # noqa: E1120 -- missing `alphabet`
    with pytest.raises(ValueError):
        d.decode(np.array([0, 1, 2]), "ACDEF")


def test_coded_bank_carries_its_alphabet_and_refuses_a_file_without_one(tmp_path):
    b = d.CodedBank.of(["ACDEF", "GHIKL"])
    p = str(tmp_path / "bank.npz")
    b.save(p)
    back = d.CodedBank.load(p)
    assert back.alphabet == d.ALPHABET
    assert back.decode_all() == ["ACDEF", "GHIKL"]
    np.savez(str(tmp_path / "bare.npz"), codes=b.codes)
    with pytest.raises(ValueError, match="without an alphabet"):
        d.CodedBank.load(str(tmp_path / "bare.npz"))


def test_encode_matches_the_repositorys_own_encoder():
    import s7.audit as audit
    for p in d.load()[:80]:
        assert (d.encode(p.seq) == audit.encode(p.seq)).all()


def test_blosum62_matches_both_hand_retyped_copies():
    """`s5/lib.py` and `s7/audit.py` each retyped BLOSUM62 by hand. If any of the three
    disagreed, every retrieval pool built from the odd one out would be different."""
    import s5.lib as lib
    import s7.audit as audit
    assert np.abs(d.BLOSUM62 - lib.B62).max() == 0.0
    assert np.abs(d.BLOSUM62 - audit.B62).max() == 0.0
    assert np.abs(d.BLOSUM62 - d.BLOSUM62.T).max() == 0.0


def test_blosum_similarity_matches_the_loop_it_replaces():
    z = np.load(_universes(1)[0], allow_pickle=True)
    S = np.asarray(z["S"])[:3000]
    q = d.encode(str(z["seq"]))
    fast = d.blosum_similarity(q, S)
    slow = np.array([sum(d.BLOSUM62[a, b] for a, b in zip(q, row)) for row in S])
    assert np.abs(fast - slow).max() < 1e-9
    # and it reproduces the deposited similarities
    assert np.abs(fast - np.asarray(z["sim"], float)[:3000]).max() < 1e-4


# 2. Identity and its normalisation
def test_identity_matches_the_legacy_alignment_exactly(peptides):
    import peptide_db as legacy
    worst = 0.0
    n = 0
    for i in range(70):
        for j in range(i + 1, 70):
            worst = max(worst, abs(d.identity(peptides[i].seq, peptides[j].seq)
                                   - legacy.identity(peptides[i].seq, peptides[j].seq)))
            n += 1
    assert n > 2000
    assert worst == 0.0


def test_identity_many_is_the_scalar_form_batched(peptides):
    """The batched DP against the per-pair one, over 31,480 real pairs. Exactly 0.0 --
    including the ragged padding, which is why the padding is right-side only."""
    import peptide_db as legacy
    seqs = [p.seq for p in peptides]
    worst = 0.0
    for q in peptides[:40]:
        a = d.identity_many(q.seq, seqs)
        b = np.array([legacy.identity(q.seq, s) for s in seqs])
        worst = max(worst, float(np.abs(a - b).max()))
    assert worst == 0.0


def test_the_batched_and_scalar_paths_agree_across_the_size_switch(peptides):
    """`_BATCH_MIN` selects an implementation for speed only, so both sides of it must give
    the same number on the same input."""
    seqs = [p.seq for p in peptides[:60]]
    q = peptides[100].seq
    small = d._nw_matches(q, seqs[:10])
    assert len(small) == 10
    big = d._nw_matches(q, seqs)
    for k in range(10):
        assert small[k] == d._nw_matches_scalar(q, seqs[k])
        assert big[k] == d._nw_matches_scalar(q, seqs[k])


def test_identity_normalises_by_the_longer_sequence_and_that_is_leaky():
    """The convention pin, and the demonstration of why it is documented as leaky.

    1CEK (13 residues) is a verbatim substring of 1A11 (25). `identity` reports 0.520,
    which passes a 0.6 filter; `containment` reports 1.000, which is the truth about the
    leak. 573 library members across 58 targets clear the shipped member-level filter while
    containing their target at >= 0.6 by this exact mechanism. The behaviour is preserved
    because results on record depend on it -- it is pinned, not endorsed.
    """
    a, b = d.by_pdb("1CEK"), d.by_pdb("1A11")
    if a is None or b is None:
        pytest.skip("1CEK / 1A11 absent from the database")
    assert a.seq in b.seq                                  # verbatim containment
    assert abs(d.identity(a.seq, b.seq) - 13.0 / 25.0) < 1e-12
    assert abs(d.identity(a.seq, b.seq) - 0.520) < 1e-9
    assert d.identity(a.seq, b.seq) < d.IDENTITY_THRESHOLD  # it PASSES the filter
    assert abs(d.containment(a.seq, b.seq) - 1.0) < 1e-12
    assert d.containment(a.seq, b.seq) >= d.IDENTITY_THRESHOLD


def test_identity_is_symmetric_and_one_on_self(peptides):
    for p in peptides[:40]:
        assert abs(d.identity(p.seq, p.seq) - 1.0) < 1e-12
        assert abs(d.containment(p.seq, p.seq) - 1.0) < 1e-12
    for a, b in zip(peptides[:30], peptides[30:60]):
        assert abs(d.identity(a.seq, b.seq) - d.identity(b.seq, a.seq)) < 1e-12


def test_the_composition_prefilter_is_admissible(peptides):
    """It must only ever OVER-estimate. The k-mer prefilter it replaced did not, and one
    benchmark target was provably trained on a 0.70-identity homolog because of it."""
    seqs = [p.seq for p in peptides]
    Cm, lens = d.composition_matrix(seqs)
    for q in peptides[:25]:
        bound = d.max_possible_identity_many(q.seq, Cm, lens)
        true = d.identity_many(q.seq, seqs)
        assert (bound + 1e-12 >= true).all(), float((true - bound).max())


def test_vectorised_prefilter_equals_the_dictionary_form(peptides):
    seqs = [p.seq for p in peptides]
    Cm, lens = d.composition_matrix(seqs)
    for q in peptides[:15]:
        fast = d.max_possible_identity_many(q.seq, Cm, lens)
        slow = np.array([d.max_possible_identity(q.seq, s) for s in seqs])
        assert np.abs(fast - slow).max() < 1e-12


# 3. DATABASE, CLUSTERS, FOLDS, TARGET SETS
def test_holdout_membership_is_identical_to_the_legacy_module(peptides):
    import peptide_db as legacy
    for p in peptides[:120]:
        assert ({x.pdb for x in d.holdout(p.seq)}
                == {x.pdb for x in legacy.holdout(p.seq)})


def test_holdout_removes_the_target_and_everything_above_the_threshold(peptides):
    for p in peptides[:30]:
        kept = d.holdout(p.seq)
        assert all(q.seq != p.seq for q in kept)
        ident = d.identity_many(p.seq, [q.seq for q in kept])
        assert (ident < d.IDENTITY_THRESHOLD).all()


def test_clusters_reproduce_the_pinned_assignment(tmp_path, monkeypatch):
    """Recompute from scratch and compare against `peptide_clusters.json`.

    The pinned file is load-bearing: fold indices are a shuffle of cluster ids, so a
    renumbering silently reassigns every sequence and a target can be scored by a model
    that trained on it, with nothing raised anywhere.
    """
    pinned = json.load(open(os.path.join(ROOT, "peptide_clusters.json")))
    monkeypatch.setattr(d, "_CLUSTER_CACHE", str(tmp_path / "cl.json"))
    fresh = d.clusters()
    assert fresh == pinned["assign"]
    assert len(set(fresh.values())) == 470


def test_folds_are_pinned_and_cover_every_sequence(peptides):
    f = d.folds(5)
    assert all(p.seq in f for p in peptides)
    assert sorted(set(f.values())) == [0, 1, 2, 3, 4]
    # a fold must be a union of whole identity clusters, or holdout leaks across folds
    cl = d.clusters()
    by_cluster = {}
    for s, c in cl.items():
        by_cluster.setdefault(c, set()).add(f[s])
    assert all(len(v) == 1 for v in by_cluster.values())


def test_benchmark_membership_matches_the_manifest():
    man = os.path.join(ROOT, "results", "benchmark_manifest.json")
    if not os.path.exists(man):
        pytest.skip("benchmark manifest absent")
    want = [t["pdb"] for t in json.load(open(man))["targets"]]
    assert [p.pdb for p in d.benchmark()] == want


def test_dev_set_is_cluster_disjoint_from_the_benchmark():
    """Sequence-disjointness is not enough: four dev peptides once shared a cluster with
    six reported targets, and hyperparameters fitted there transferred to the homologs and
    not the rest -- a difference-in-differences of 0.424 A."""
    cl = d.clusters()
    bench = {cl[p.seq] for p in d.benchmark()}
    dev = {cl[p.seq] for p in d.dev_set()}
    assert not (bench & dev)


def test_target_sets_match_the_legacy_module():
    import peptide_db as legacy
    assert [p.pdb for p in d.benchmark()] == [p.pdb for p in legacy.benchmark()]
    assert [p.pdb for p in d.dev_set()] == [p.pdb for p in legacy.dev_set()]
    assert [p.pdb for p in d.benchmark_set()] == [p.pdb for p in legacy.benchmark_set()]
    if os.path.exists(d.MONOMER_MANIFEST):
        assert ([p.pdb for p in d.monomer_benchmark()]
                == [p.pdb for p in legacy.monomer_benchmark()])


# 4. RETRIEVAL TIE-BREAKING
def test_stable_argsort_reproduces_every_deposited_pool():
    """The tie-break pin. `s7/audit.py` used a stable argsort and `s5/inband.py` did not.

    Over all 126 universes the stable order reproduces the deposited ``order`` array with
    zero mismatches. The unstable sort differs on real targets by as many as 47 of 500
    members -- a different pool, and therefore a different result. Leaving that to the
    sorting algorithm is not a defensible choice, so it is fixed.
    """
    worst_stable = 0
    worst_unstable = 0
    for f in _universes():
        z = np.load(f, allow_pickle=True)
        sim = np.asarray(z["sim"], float)
        want = np.asarray(z["order"])[:500]
        worst_stable = max(worst_stable, int((d.top_k(sim, 500) != want).sum()))
        unstable = np.argsort(-sim)[:500]
        worst_unstable = max(
            worst_unstable,
            len(set(unstable.tolist()) ^ set(want.tolist())) // 2)
    assert worst_stable == 0
    assert worst_unstable > 0, "the pin would be vacuous if ties never mattered"


def test_windows_are_contiguous_slices_of_their_parents(peptides):
    pool = [p for p in peptides[:60] if p.n >= 12]
    W, S, src = d.windows(pool, 9)
    for k in (0, 5, len(src) // 2, len(src) - 1):
        m, s = src[k]
        assert np.abs(W[k] - pool[m].ca[s:s + 9]).max() == 0.0
        assert d.decode(S[k], d.ALPHABET) == pool[m].seq[s:s + 9]


def test_windows_carry_their_parents_torsions(peptides):
    pool = [p for p in peptides[:40] if p.n >= 12]
    W, S, src, PHI, PSI = d.windows(pool, 9, with_torsions=True)
    m, s = src[7]
    assert np.abs(PHI[7] - pool[m].phi[s:s + 9].astype(np.float32)).max() == 0.0
    assert np.abs(PSI[7] - pool[m].psi[s:s + 9].astype(np.float32)).max() == 0.0


# 5. REPRESENTATIONS
def test_torsion_representation_decodes_identically_to_the_legacy_module(peptides):
    import representations as legacy
    rng = np.random.default_rng(0)
    for p in peptides[:40]:
        for ns in (4, 8):
            a = d.make_representation("torsion", p.n, ns, sequence=p.seq)
            b = legacy.make_representation("torsion", p.n, ns, sequence=p.seq)
            assert a.n_bits == b.n_bits and a.classes == b.classes
            for _ in range(4):
                bs = b.random_bitstring(rng)
                pa, qa = a.decode(bs)
                pb, qb = b.decode(bs)
                assert np.abs(pa - pb).max() == 0.0
                assert np.abs(qa - qb).max() == 0.0
            assert a.native_bitstring(p.phi, p.psi) == b.native_bitstring(p.phi, p.psi)


def test_legacy_library_encoding_is_preserved_bit_for_bit(peptides):
    """`run_8state_seed0.py` and `amber_obc2.py` are pinned to persisted bitstring->energy
    pairs verified to max absolute error 0.0. Changing the library changes the coordinates
    a bitstring decodes to and therefore every energy in that set."""
    import representations as legacy
    rng = np.random.default_rng(1)
    for p in peptides[:20]:
        for ns in (4, 8):
            a = d.make_representation("torsion", p.n, ns, sequence=p.seq,
                                      legacy_library=True)
            b = legacy.make_representation("torsion", p.n, ns, sequence=p.seq,
                                           legacy_library=True)
            assert a.n_bits == b.n_bits and a.n_chi_bits == 0
            bs = b.random_bitstring(rng)
            assert np.abs(np.asarray(a.decode(bs)) - np.asarray(b.decode(bs))).max() == 0.0


def test_lattice_decode_is_identical_and_in_angstroms(peptides):
    import representations as legacy
    rng = np.random.default_rng(2)
    for p in peptides[:30]:
        a = d.make_representation("lattice", p.n, sequence=p.seq)
        b = legacy.make_representation("lattice", p.n, sequence=p.seq)
        for _ in range(5):
            bs = b.random_bitstring(rng)
            assert np.abs(a.decode(bs) - b.decode(bs)).max() == 0.0
        ca = a.decode(b.random_bitstring(rng))
        step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        assert np.abs(step - d.CA_CA_DISTANCE).max() < 1e-9


def test_state_zero_is_helical_and_state_one_extended_in_every_library():
    """A convention callers rely on to build reference structures: `[0] * n` must be "as
    helical as this sequence can be" whatever the composition."""
    for ns, libs in d.STATE_LIBRARIES.items():
        for name, lib in libs.items():
            assert lib[0] == (-63.0, -42.0) or lib[0][0] < 0 and lib[0][1] < 0, (ns, name)
            assert lib[1][1] > 90.0, (ns, name)


def test_batched_decode_matches_the_per_bitstring_form(peptides):
    rng = np.random.default_rng(4)
    p = peptides[0]
    rep = d.make_representation("torsion", p.n, 4, sequence=p.seq)
    states = rng.integers(0, 4, size=(32, p.n))
    phi, psi = rep.decode_batch(states)
    for k in range(32):
        bs = rep.bitstring_from_states(states[k])
        a, b = rep.decode(bs)
        assert np.abs(phi[k] - a).max() == 0.0
        assert np.abs(psi[k] - b).max() == 0.0


# 6. THE CACHE
def test_key_changes_with_every_parameter_and_with_its_type():
    k = C.key("ns", 1, a=1, b="x")
    assert C.key("ns", 1, a=1, b="x") == k                 # deterministic
    assert C.key("ns", 1, a=2, b="x") != k                 # value
    assert C.key("ns", 2, a=1, b="x") != k                 # version
    assert C.key("other", 1, a=1, b="x") != k              # namespace
    assert C.key("ns", 1, a=1) != k                        # a MISSING parameter
    assert C.key("ns", 1, a=1.0, b="x") != k               # int vs float
    assert C.key("ns", 1, a=True, b="x") != k              # bool vs int
    assert C.key("ns", 1, a="1", b="x") != k               # str vs int
    assert (C.key("ns", 1, a=np.zeros(3, np.float32))
            != C.key("ns", 1, a=np.zeros(3, np.float64)))  # dtype
    assert (C.key("ns", 1, a=(1, 2)) != C.key("ns", 1, a=[1, 2]))   # tuple vs list


def test_unhashable_parameters_are_refused_rather_than_stringified():
    """`str()` of a plain object embeds its address, so a stringified key would be
    nondeterministic across runs -- and could collide on an address reuse."""
    class Thing:
        pass
    with pytest.raises(TypeError):
        C.key("ns", 1, a=Thing())


def test_store_and_load_round_trip(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "ROOT", str(tmp_path))
    C._MEM.clear()
    del C._MEM_ORDER[:]
    arr = np.arange(12, dtype=np.float64).reshape(3, 4)
    C.store("t", 1, {"x": arr}, seq="ABC", k=3)
    C._MEM.clear()
    del C._MEM_ORDER[:]
    got = C.load("t", 1, seq="ABC", k=3)
    assert got is not None and np.abs(got["x"] - arr).max() == 0.0
    assert C.load("t", 1, seq="ABD", k=3) is None
    assert C.load("t", 2, seq="ABC", k=3) is None


def test_a_key_that_omits_a_parameter_is_caught_not_served(tmp_path, monkeypatch):
    """The whole point. A digest collision is astronomically unlikely; a bug in the key
    construction is not, and the sidecar is what turns it from silent into loud."""
    monkeypatch.setattr(C, "ROOT", str(tmp_path))
    C._MEM.clear()
    del C._MEM_ORDER[:]
    C.store("t", 1, {"x": np.ones(2)}, seq="ABC")
    k = C.key("t", 1, seq="ABC")
    npz, side = C._paths("t", k)
    meta = json.load(open(side))
    meta["params"] = json.dumps(C._canon({"seq": "ABC", "forgotten": 7}),
                                sort_keys=True, separators=(",", ":"))
    json.dump(meta, open(side, "w"))
    C._MEM.clear()
    del C._MEM_ORDER[:]
    with pytest.raises(C.CacheCollision):
        C.load("t", 1, seq="ABC")


def test_cached_computes_once_then_serves(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "ROOT", str(tmp_path))
    C._MEM.clear()
    del C._MEM_ORDER[:]
    calls = []

    def compute():
        calls.append(1)
        return {"v": np.array([42.0])}

    for _ in range(3):
        assert C.cached("t", 1, compute, q=1)["v"][0] == 42.0
    assert len(calls) == 1


def test_cache_can_be_disabled(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "ROOT", str(tmp_path))
    monkeypatch.setattr(C, "DISABLED", True)
    C._MEM.clear()
    del C._MEM_ORDER[:]
    C.store("t", 1, {"x": np.ones(2)}, seq="Z")
    C._MEM.clear()
    del C._MEM_ORDER[:]
    assert C.load("t", 1, seq="Z") is None


def test_fold_fragments_matches_the_legacy_filter():
    """Cached, but it must be the same set: it decides what a fold model may train on."""
    import distogram as legacy
    a = {p.pdb for p in d.fold_fragments(0, 5)}
    b = {p.pdb for p in legacy._fold_fragments(0, 5)}
    assert a == b


def test_prediction_features_still_agree_with_priors(peptides):
    """`core.predict` owns the pair-feature stack that `priors.py` used to.

    `priors.py` is out of the consolidated import closure -- `core.predict` replaces
    `distogram`, which was its only production consumer -- but it stays on disk for
    `DistogramPrior` and `TorsionMRF`, which are a different scoring form. So the property
    table and the pair features exist twice, and the ONLY thing stopping them drifting is
    this test. Difference is exactly 0.0 over 400 real sequences.
    """
    from core import predict as cp
    import priors
    assert np.abs(cp.BIN_EDGES - priors.BIN_EDGES).max() == 0.0
    assert cp.NBINS == priors.NBINS
    worst = 0.0
    for p in peptides[:400]:
        assert np.abs(cp.residue_props(p.seq) - priors._prop(p.seq)).max() == 0.0
        a = priors.pair_features(p.seq)
        b = cp.pair_features(p.seq)
        worst = max(worst, float(np.abs(a[0] - b[0]).max()))
        assert (a[1] == b[1]).all() and (a[2] == b[2]).all()
    assert worst == 0.0, worst


def test_distogram_predictions_are_identical_end_to_end():
    """The whole prediction stack, against `distogram.py`, on trained checkpoints.

    Restricted to sequences already in the 21 MB hot ESM cache so that running the test
    suite cannot itself pull the 1.5 GB bank into a long-lived process -- which is exactly
    the failure this consolidation exists to remove, and which was reproduced today by
    calling the legacy `pairnet.residue_features` on a cold sequence.
    """
    if not os.path.exists(os.path.join(ROOT, "esm_small.npz")):
        pytest.skip("hot ESM cache absent")
    if not glob.glob(os.path.join(ROOT, "distogram_models", "fold*_esm_frag.pt")):
        pytest.skip("no trained distogram checkpoints")
    from core import predict as cp
    import distogram as legacy
    hot = {str(s) for s in np.load(os.path.join(ROOT, "esm_small.npz"),
                                   allow_pickle=True)["seqs"]}
    targets = [p for p in d.load() if p.seq in hot][:4]
    assert targets
    for p in targets:
        a = legacy.Distogram.for_target(p.seq)
        b = cp.Distogram.for_target(p.seq)
        assert np.abs(a.prob - b.prob).max() == 0.0
        assert np.abs(a.expected - b.expected).max() == 0.0
        W = np.stack([q.ca for q in d.load() if q.n == p.n][:32])
        if len(W) > 2:
            assert np.abs(np.asarray(a.score(W)) - np.asarray(b.score(W))).max() == 0.0


def test_feature_width_never_needs_the_big_esm_bank():
    """It is an integer. Learning it must not cost 1.5 GB of resident memory."""
    from core import predict as cp
    w = cp.feature_width(True)
    assert w in (42, 183), w
    assert cp.feature_width(True) == w          # and it is served, not recomputed
    assert cp.feature_width(False) == 42


def test_memory_gate_reads_in_process():
    """`Get-CimInstance` costs minutes per call on a loaded box, which turns the gate that
    protects the box into the slowest thing on it."""
    m = C.mem_pct()
    assert m is None or 0 <= m <= 100
