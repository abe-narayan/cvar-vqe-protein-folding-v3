"""s25/resultslab/test_export.py -- the export layer's unit tests.

Run:  python -m s25.resultslab.test_export        (no pytest dependency)
      pytest s25/resultslab/test_export.py        (also works)

The load-bearing test is `test_roundtrip_against_instrument`: an exported structure must
reproduce, when read back off disk and scored against the exported native, the RMSD that
`s12.instrument.ca_rmsd` computes in memory.  Everything else in the results lab is
decoration if that one fails.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                              # noqa: E402
from s25.resultslab import exportlib as EX                   # noqa: E402

#: A spread of the pinned dev set: shortest, longest, a FAIL18 member, two mid targets.
SAMPLE = ["1A13", "1A1P", "1CB3", "1ID6", "2BFI", "9WXK"]


# ------------------------------------------------------------------ the T-number mapping
def test_target_map_is_stable_and_total():
    m = EX.build_target_map()
    assert m["n"] == 126, m["n"]
    assert m["complete"] is True
    pdbs = [t["pdb"] for t in I.targets()]
    assert [r["pdb_id"] for r in m["targets"]] == pdbs
    assert m["t_of_pdb"][pdbs[0]] == "T001"
    assert m["pdb_of_t"]["T126"] == pdbs[125]
    assert len(set(m["t_of_pdb"].values())) == 126
    # idempotent: a second build must return the same object, not renumber
    again = EX.build_target_map()
    assert again["digest"] == m["digest"]


def test_stem_accepts_both_id_forms_and_rejects_unknown_config():
    assert EX.stem("1A13", "legacy") == "T001__legacy"
    assert EX.stem("T001", "legacy") == "T001__legacy"
    assert EX.stem("1A13", "native") == "T001__native"
    assert EX.stem("1A13", "_noise") == "T001___noise"
    try:
        EX.stem("1A13", "not_a_configuration")
    except KeyError:
        pass
    else:
        raise AssertionError("unknown configuration was accepted")


def test_provenance_travels_with_the_file():
    """GATE (c), the primary one.  A structure must declare itself IN ITS OWN HEADER.

    The quarantined proof build wrote 1,134 PDBs of native-plus-noise, six of seven
    configurations, and not one said so in its header -- `is_synthetic` reached the JSON and
    never the artefact.  Read on its own, `T001__distogram.pdb` announced RMSD_CA 0.000000
    as a genuine result.
    """
    rec = I.shipped_record("1A13")
    pc = np.asarray(rec["avg_ca"], float)

    good = EX.export_prediction("1A13", "_incumbent", pc, basis=EX.SECONDARY_BASIS,
                                provenance=EX.PROV_GENUINE)
    hdr = EX.verify_header(good["path"])
    for k in EX.REQUIRED_HEADER_KEYS:
        assert str(hdr.get(k, "")).strip(), k
    assert hdr["provenance"] == EX.PROV_GENUINE
    assert len(hdr["module_sha"]) >= 8 and len(hdr["git_commit"]) >= 7

    # a synthetic file must SAY SO, and the frozen build must refuse it
    synth = EX.export_prediction("1A13", "_noise", pc, basis=EX.SECONDARY_BASIS,
                                 provenance=EX.PROV_SYNTHETIC)
    txt = open(synth["path"]).read()
    assert "REMARK 999 PROVENANCE synthetic" in txt
    assert "NOT A RESULT -- TEST FIXTURE, SYNTHETIC" in txt
    for allowed in (EX.PROV_GENUINE,):
        try:
            EX.verify_header(synth["path"], allow_provenance=(allowed,))
        except EX.ProvenanceViolation as e:
            assert "PROVENANCE GATE FAILED" in str(e)
        else:
            raise AssertionError("a synthetic file passed the frozen build's header gate")

    # undeclared provenance is UNKNOWN and also fails
    unk = EX.export_prediction("1A13", "_zero", pc, basis=EX.SECONDARY_BASIS)
    assert EX.read_pdb(unk["path"])["meta"]["provenance"] == EX.PROV_UNKNOWN
    try:
        EX.verify_header(unk["path"])
    except EX.ProvenanceViolation:
        pass
    else:
        raise AssertionError("an undeclared file passed the header gate")

    # a header key stripped out of the file must be detected, not assumed
    stripped = str(synth["path"]) + ".stripped.pdb"
    with open(stripped, "w", newline="\n") as fh:
        fh.writelines(ln for ln in open(synth["path"]) if "REMARK 999 PROVENANCE" not in ln)
    try:
        EX.verify_header(stripped, allow_provenance=(EX.PROV_GENUINE, EX.PROV_SYNTHETIC))
    except EX.ProvenanceViolation as e:
        assert "carries no provenance" in str(e)
    else:
        raise AssertionError("a header with PROVENANCE removed passed the gate")
    os.remove(stripped)


def test_export_is_atomic_and_never_leaves_pending():
    """One write, tmp + os.replace.  `T006__amber.pdb` was left with `RMSD_CA PENDING`."""
    import glob
    rec = I.shipped_record("1A1P")
    r = EX.export_prediction("1A1P", "_incumbent", np.asarray(rec["avg_ca"], float),
                             basis=EX.SECONDARY_BASIS, provenance=EX.PROV_GENUINE)
    txt = open(r["path"]).read()
    assert "PENDING" not in txt
    assert not glob.glob(os.path.join(os.path.dirname(r["path"]), "*.tmp*"))
    # the header RMSD is the RMSD of the bytes, computed BEFORE the single write
    got = EX.read_pdb(r["path"])
    # the header carries %.6f, so agreement is exact to half a unit in the last place
    assert abs(I.ca_rmsd(got["CA"], EX.native_ca("1A1P")) - float(got["meta"]["rmsd_ca"])) < 5e-7


def test_quantise_matches_the_bytes_exactly():
    """`quantise` must agree with the file bit for bit -- np.round would not, on ties."""
    rng = np.random.default_rng(3)
    seq = EX.sequence_of("1CB3")
    coords = {"CA": rng.normal(0, 5, (len(seq), 3))}
    q = EX.quantise(seq, coords)
    p_ = os.path.join(EX.FIXTURE_DIR, "_tmp_q.pdb")
    EX.write_structure(p_, seq, coords, {"basis": EX.SECONDARY_BASIS})
    got = EX.read_pdb(p_)
    os.remove(p_)
    assert np.array_equal(q["CA"], got["CA"])


def test_difficulty_gate_separates_genuine_from_noise():
    """GATE (b), the one that actually discriminates.  Calibration re-measured here."""
    pb = EX.pool_best()
    names = [t["pdb_id"] for t in EX.target_rows()]
    inc = {p: I.ca_rmsd(np.asarray(I.shipped_record(p)["ca"], float), EX.native_ca(p))
           for p in names}
    g = EX.difficulty_gate(inc, "production")
    assert g["status"] == "PASS", g
    assert g["corr_with_pool_best"] > 0.6 and g["sd"] > 1.5, g

    import zlib
    for sigma in (0.9, 1.6, 2.4):
        syn = {}
        for p in names:
            nat = EX.native_ca(p)
            rng = np.random.default_rng(zlib.crc32(("%s|g|0" % p).encode()) & 0xFFFFFFFF)
            syn[p] = I.ca_rmsd(nat + rng.normal(0, sigma, nat.shape), nat)
        soft = EX.difficulty_gate(syn, "legacy", hard=False)
        assert soft["status"] == "FAIL", (sigma, soft)
        assert abs(soft["corr_with_pool_best"]) < 0.25, soft
        try:
            EX.difficulty_gate(syn, "legacy")
        except EX.DifficultyGateViolation as e:
            assert "DIFFICULTY GATE FAILED" in str(e)
        else:
            raise AssertionError("the difficulty gate passed isotropic noise at sigma=%s" % sigma)
    return g


def test_production_arm_is_lam_0_3_and_reproduces_the_shipped_synthesis():
    """We ship SYNTHESIS stage 3 (`ca`, 3.2148 A), not the lam=0 fit arm (`fit_ca`, 3.2041)."""
    assert EX.PRODUCTION_LAM == 0.3
    for pdb in ["1A13", "1A1P", "1CB3"]:
        rec = I.shipped_record(pdb)
        ch = EX.chain_from_ca(np.asarray(rec["avg_ca"], float), EX.sequence_of(pdb),
                              EX.fold_of(pdb))                      # default = PRODUCTION_LAM
        assert float(np.abs(ch["ca"] - np.asarray(rec["ca"], float)).max()) == 0.0
        assert abs(I.ca_rmsd(ch["ca"], EX.native_ca(pdb)) - rec["rmsd_arm"]) < 1e-9


def test_fixtures_are_quarantined_and_cannot_look_like_results():
    """Fixture structures may never land beside real configuration exports.

    The round-trip fixture used to name its controls "distogram" (the native, RMSD exactly
    0.000000) and "legacy" (native + noise), so the test's own output contained the line
    `T001  distogram  14  0.000000  0.000000`.  One careless paste of that is the most
    damaging sentence this project could emit.  Controls now carry `_`-prefixed names and
    land in their own directory.
    """
    root = os.path.abspath(EX.RESULTS_DIR)
    for f in EX.FIXTURES:
        assert f.startswith("_") and f not in EX.CONFIGURATIONS
        p_ = os.path.abspath(EX.struct_path("1A13", f, EX.PRIMARY_BASIS))
        # GATE (d): fixtures live in the SYSTEM TEMP DIRECTORY, never under results/
        assert not p_.startswith(root), p_
        assert os.path.abspath(EX.FIXTURE_DIR).startswith(
            os.path.abspath(__import__("tempfile").gettempdir()))
    real = os.path.abspath(EX.struct_path("1A13", "legacy", EX.PRIMARY_BASIS))
    assert real.startswith(root) and "_fixtures" not in real.replace("\\", "/")


def test_synthetic_provider_cannot_take_a_real_configuration_name():
    from s25.resultslab import providers as PV
    PV.clear()
    try:
        PV.register("legacy", PV.synthetic(sigma=1.0, seed=0, tag="x"))
    except ValueError as e:
        assert "REFUSED" in str(e)
    else:
        PV.clear()
        raise AssertionError("a synthetic provider was bound to a real configuration name")
    PV.register("_noise", PV.synthetic(sigma=1.0, seed=0, tag="x"))   # fixture label: allowed
    PV.clear()


def test_absent_provenance_is_unknown_not_genuine(tmpdir=None):
    """`from_json` must not stamp a file with no provenance as REAL."""
    import json as _json
    from s25.resultslab import providers as PV
    nat = EX.native_ca("1A13").tolist()
    p_ = os.path.join(EX.SUMMARY_DIR, "_tmp_prov.json")
    with open(p_, "w") as fh:
        _json.dump({"1A13": nat}, fh)
    got = PV.from_json(p_)("1A13")
    os.remove(p_)
    assert got["is_synthetic"] == PV.UNKNOWN, got["is_synthetic"]
    assert got["is_synthetic"] is not False


def test_pool_oracle_gate_fires_on_a_leaking_configuration():
    """THE RELEASE GATE.  A configuration that beats its own candidate pool must not ship."""
    pb = EX.pool_best()
    sample = {k: pb[k] for k in list(sorted(pb))[:40]}
    clean = {k: v + 1.5 for k, v in sample.items()}
    ok = EX.pool_oracle_gate(clean, "legacy")
    assert ok["status"] == "PASS" and not ok["aggregate_violation"]
    leaking = {k: v * 0.25 for k, v in sample.items()}          # e.g. a provider on nat_ca
    try:
        EX.pool_oracle_gate(leaking, "legacy")
    except EX.PoolOracleViolation as e:
        assert "POOL-ORACLE GATE FAILED" in str(e)
    else:
        raise AssertionError("the pool-oracle gate did not fire on a leaking configuration")
    soft = EX.pool_oracle_gate(leaking, "legacy", hard=False)
    assert soft["status"] == "FAIL" and soft["n_per_target_violations"] > 0
    # the pinned mean is what the gate is calibrated against
    assert abs(float(np.mean(list(pb.values()))) - EX.POOL_BEST_MEAN) < 2e-3


def test_point_cloud_is_not_a_molecule():
    """The measurement behind the basis ruling, RECOMPUTED here rather than quoted.

    Over all 126 dev targets this module measures, independently of the ledger:
        point cloud  mean virtual Ca-Ca bond 2.9614 A, shortest single bond 0.6492 A,
                     80/126 targets under 3.4 A  -> 22.3% contracted
        built chain  3.8040 A        native  3.8122 A
    The full-set figures are printed by `--full`; this fast test asserts the DIRECTION on
    the 6-target sample, whose mean (3.5059 A) is milder than the 126-target mean because
    the sample happens to miss the worst offenders.  Asserting the 126-target number on 6
    targets would be a test that passes for the wrong reason.
    """
    pc, bc, nb, mn = [], [], [], []
    for pdb in SAMPLE:
        rec = I.shipped_record(pdb)
        a = EX.ca_bond_stats(np.asarray(rec["avg_ca"], float))
        pc.append(a["ca_bond_mean"]); mn.append(a["ca_bond_min"])
        bc.append(EX.ca_bond_stats(np.asarray(rec["fit_ca"], float))["ca_bond_mean"])
        nb.append(EX.ca_bond_stats(EX.native_ca(pdb))["ca_bond_mean"])
    mpc, mbc, mnb = map(np.mean, (pc, bc, nb))
    assert mpc < mnb - 0.15, (mpc, mnb)     # the point cloud is contracted, on every sample
    assert min(mn) < 3.0                    # and carries bonds no molecule could have
    assert abs(mbc - mnb) < 0.05, (mbc, mnb)   # the built chain has native-like virtual bonds


# ------------------------------------------------------------------ writer equivalence
def test_atom_lines_match_core_writer():
    """`_atom_lines` must be byte-identical to `core.geometry.write_pdb`'s body.

    This is why the lab does not carry a second, divergent PDB format.
    """
    from core import geometry as geo
    rng = np.random.default_rng(0)
    seq = EX.sequence_of("1A13")
    n = len(seq)
    coords = {k: rng.normal(0, 5, (n, 3)) for k in ("N", "CA", "C", "O")}
    ref_path = os.path.join(EX.SUMMARY_DIR, "_tmp_ref.pdb")
    geo.write_pdb(ref_path, seq, coords, remark="x")
    with open(ref_path) as fh:
        ref = [ln for ln in fh if ln.startswith("ATOM  ")]
    os.remove(ref_path)
    mine = EX._atom_lines(seq, coords)
    assert mine == ref, "ATOM formatting diverged from core.geometry.write_pdb"

    # CA-only export really is CA-only
    ca_only = EX._atom_lines(seq, {"CA": coords["CA"]})
    assert len(ca_only) == n
    assert all(ln[12:16] == " CA " for ln in ca_only)


def test_reader_matches_core_parser_on_a_backbone_file():
    """The local column reader must agree with `core.geometry.parse_pdb` where both apply."""
    from core import geometry as geo
    rng = np.random.default_rng(1)
    seq = EX.sequence_of("1A1P")
    n = len(seq)
    coords = {k: rng.normal(0, 5, (n, 3)) for k in ("N", "CA", "C", "O")}
    p = os.path.join(EX.SUMMARY_DIR, "_tmp_bb.pdb")
    EX.write_structure(p, seq, coords, {"basis": EX.BASIS_BUILT_CHAIN})
    mine = EX.read_pdb(p)
    got_seq, N, CA, C = geo.parse_pdb_ensemble(p, use_cache=False)[0]
    os.remove(p)
    assert mine["sequence"] == got_seq == seq
    # core.geometry round-trips deposited coordinates through float32 by design; compare at
    # the quantisation scale, not bitwise.
    assert np.abs(mine["CA"] - CA).max() < 1e-3
    assert np.abs(mine["N"] - N).max() < 1e-3


# ------------------------------------------------------------------ THE round-trip
def _roundtrip_rows(pdbs=SAMPLE, seed=7):
    """Three FIXTURE arms per target: the incumbent replay, the native itself, and noise.

    The labels are `_incumbent`, `_zero`, `_noise` and NOT real configuration names.  When
    they were `production`/`distogram`/`legacy`, this function's own printed table contained
    `T001  distogram  14  0.000000  0.000000` -- a line that reads as a perfect result for a
    real method.  The arms are unchanged; only the labels are safe now.
    """
    rows = []
    rng = np.random.default_rng(seed)
    for pdb in pdbs:
        nat = EX.native_ca(pdb)
        rec = I.shipped_record(pdb)
        arms = {
            "_incumbent": np.asarray(rec["avg_ca"], float),      # the shipped point cloud
            "_zero": np.asarray(nat, float).copy(),              # a zero-RMSD arm
            "_noise": np.asarray(nat, float) + rng.normal(0, 1.5, nat.shape),
        }
        for cfg, ca in arms.items():
            rows.append(EX.verify_roundtrip(pdb, cfg, ca, basis=EX.SECONDARY_BASIS))
    return rows


def test_roundtrip_against_instrument():
    rows = _roundtrip_rows()
    worst = max(r["abs_err_vs_instrument"] for r in rows)
    worst_c = max(r["max_coord_abs_err"] for r in rows)
    assert all(r["ok"] for r in rows), [r for r in rows if not r["ok"]]
    assert worst < EX.RMSD_TOL, worst
    assert worst_c <= 5.001e-4, worst_c
    # the exact-native arm must come back as a genuine zero, not merely small
    zeros = [r for r in rows if r["configuration"] == "_zero"]
    assert max(r["rmsd_from_file_vs_file_native"] for r in zeros) < 1e-6
    return rows


def test_roundtrip_on_the_built_chain_basis():
    """The PRODUCTION basis, where the file holds N/CA/C/O.  Full-instrument version is
    `verify_all_built_chain`, run by `--full`; this is the fast gate."""
    rows = []
    for pdb in ["1A13", "1A1P", "1CB3"]:
        rec = I.shipped_record(pdb)
        ch = EX.chain_from_ca(np.asarray(rec["avg_ca"], float),
                              EX.sequence_of(pdb), EX.fold_of(pdb), lam=0.0)
        r = EX.verify_roundtrip(pdb, "_incumbent", ch["ca"],
                                basis=EX.BASIS_BUILT_CHAIN, backbone=ch["backbone"])
        # the lam=0 fit arm must reproduce the cached `fit_ca` this repo already shipped
        r["vs_cached_fit_ca"] = float(np.abs(ch["ca"] - np.asarray(rec["fit_ca"], float)).max())
        assert r["vs_cached_fit_ca"] == 0.0, r["vs_cached_fit_ca"]
        rows.append(r)
    assert all(r["ok"] for r in rows), rows
    return rows


def test_export_pair_emits_both_bases_from_one_call():
    """`export_pair` is how "always paired" becomes a property of the code."""
    rec = I.shipped_record("1A13")
    pair = EX.export_pair("1A13", "_incumbent", np.asarray(rec["avg_ca"], float))
    assert pair["basis"] == EX.BASIS_BUILT_CHAIN
    assert pair["basis_secondary"] == EX.BASIS_POINT_CLOUD
    assert abs(pair["rmsd"] - rec["rmsd_arm"]) < EX.RMSD_TOL     # SYNTHESIS, lam=0.3
    assert abs(pair["rmsd_secondary"] - rec["rmsd_avg"]) < EX.RMSD_TOL
    assert pair["delta"] > 0                       # the honest number is the worse one
    for which in ("primary", "secondary"):
        got = EX.read_pdb(pair[which]["path"])
        assert got["meta"]["basis"] == pair[which]["basis"]
        assert abs(float(got["meta"]["rmsd_ca"]) - pair[which]["rmsd"]) < 1e-6
    return pair


def test_absent_is_absent_not_estimated():
    """A missing prediction must raise, never be silently filled."""
    try:
        EX.export_prediction("1A13", "_noise", np.zeros((3, 3)),
                             basis=EX.SECONDARY_BASIS)
    except ValueError as e:
        assert "shape" in str(e)
    else:
        raise AssertionError("a wrong-length prediction was accepted")
    try:
        EX.export_prediction("1A13", "_noise",
                             np.full((EX.length_of("1A13"), 3), np.nan),
                             basis=EX.SECONDARY_BASIS)
    except ValueError as e:
        assert "finite" in str(e)
    else:
        raise AssertionError("a non-finite prediction was accepted")


# ------------------------------------------------------------------ basis discipline
def test_synthetic_provider_is_deterministic():
    """Deterministic result loading is a hard constraint; a hash-seeded RNG is not.

    Python salts string hashing per process, so `hash()` in a seed made a "reproducible"
    build emit different coordinates on every run.  This asserts a fresh interpreter with a
    DIFFERENT PYTHONHASHSEED reproduces the same structure bit for bit.
    """
    import subprocess
    from s25.resultslab import providers as PV
    a = PV.synthetic(sigma=1.5, seed=3, tag="legacy")("1A13")["ca"]
    code = ("import sys; sys.path.insert(0, " + repr(ROOT) + ");"
            "from s25.resultslab import providers as PV;"
            "print(format(float(PV.synthetic(sigma=1.5, seed=3, tag='legacy')"
            "('1A13')['ca'].sum()), '.17g'))")
    outs = set()
    for hs in ("0", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=hs)
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                           env=env, cwd=ROOT)
        assert r.returncode == 0, r.stderr
        outs.add(r.stdout.strip())
    assert len(outs) == 1, "synthetic provider is not hash-seed stable: %s" % outs
    assert abs(float(next(iter(outs))) - float(a.sum())) < 1e-9


def test_basis_constants_are_distinct_and_recorded():
    assert EX.BASIS_POINT_CLOUD != EX.BASIS_BUILT_CHAIN != EX.BASIS_NATIVE
    assert EX.PRIMARY_BASIS == EX.BASIS_BUILT_CHAIN     # s25 L4/L8: the chain is production
    assert EX.SECONDARY_BASIS == EX.BASIS_POINT_CLOUD
    rec = I.shipped_record("1A13")
    r = EX.export_prediction("1A13", "_incumbent", np.asarray(rec["avg_ca"], float),
                             basis=EX.SECONDARY_BASIS)
    got = EX.read_pdb(r["path"])
    assert got["meta"]["basis"] == EX.BASIS_POINT_CLOUD
    assert "NON-PHYSICAL" in got["meta"]["basis_role"]
    assert abs(float(got["meta"]["rmsd_ca"]) - r["rmsd"]) < 1e-6


def full_built_chain_report():
    """`--full`: the round trip on ALL 126 targets, on the PRODUCTION basis."""
    print("\n  FULL-INSTRUMENT ROUND TRIP, built-chain basis, 126 targets ...")
    out = EX.verify_all_built_chain()
    print("    all_ok                              %s" % out["all_ok"])
    print("    worst |file - instrument|  chain    %.3e A  (tol %.1e)"
          % (out["worst_abs_err"], out["tolerance"]))
    print("    worst |file - instrument|  cloud    %.3e A" % out["worst_abs_err_point_cloud"])
    print("    reprojection vs the cached arm      %.3e A  (0 = bit-exact stage 3b)"
          % out["worst_reprojection_vs_cached_arm"])
    print("    mean RMSD  BUILT CHAIN (production) %.4f A" % out["mean_rmsd_built_chain"])
    print("    mean RMSD  point cloud (non-phys.)  %.4f A" % out["mean_rmsd_point_cloud"])
    print("    basis gap (chain - cloud)           %+.4f A" % out["basis_gap"])
    print("    mean virtual CA-CA bond   native    %.4f A" % out["mean_ca_bond_native"])
    print("                              chain     %.4f A" % out["mean_ca_bond_built_chain"])
    print("                              cloud     %.4f A  <- 22%% contracted, not a molecule"
          % out["mean_ca_bond_point_cloud"])
    print("    shortest virtual bond   point cloud  %.4f A  <- shorter than a covalent C-C"
          % out["min_ca_bond_point_cloud"])
    print("                            built chain  %.4f A" % out["min_ca_bond_built_chain"])
    assert out["all_ok"], "full-instrument round trip FAILED"
    return out


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    rows = None
    for t in tests:
        got = t()
        if t.__name__ == "test_roundtrip_against_instrument":
            rows = got
        print("  PASS  %s" % t.__name__)
    if rows:
        print("\n  ROUND-TRIP: exported file -> reread -> ca_rmsd  vs  instrument ca_rmsd")
        print("  %-6s %-14s %-6s %12s %12s %12s %10s"
              % ("target", "configuration", "n", "instrument", "from file", "abs err", "max dxyz"))
        for r in rows:
            print("  %-6s %-14s %-6d %12.6f %12.6f %12.2e %10.2e"
                  % (r["target_id"], r["configuration"], r["n"], r["rmsd_instrument"],
                     r["rmsd_from_file_vs_file_native"], r["abs_err_vs_instrument"],
                     r["max_coord_abs_err"]))
        print("  worst abs err %.3e A  (tolerance %.1e, quantisation bound 8.66e-04)"
              % (max(r["abs_err_vs_instrument"] for r in rows), EX.RMSD_TOL))
    print("\n  all export tests PASS")
    if "--full" in sys.argv:
        full_built_chain_report()


if __name__ == "__main__":
    main()
