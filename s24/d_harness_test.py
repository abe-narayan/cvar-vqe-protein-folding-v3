"""S24 / LANE D -- unit tests for the integration harness. Run this before either lane uses it.

    python s24/d_harness_test.py

Tests the things that would silently corrupt a lane's result rather than crash it: the two
arms seeing different inputs, the readout being applied in a different frame per arm, a
padding state reaching a readout, a basis mix surviving aggregation, an oracle label leaking
into a selection, the null not matching the operator's space, and the completion flag
passing on a row that is missing keys.
"""
from __future__ import annotations

import json
import os
import sys
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I          # noqa: E402
from s24 import d_harness as H           # noqa: E402

FAILED: list = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:
        FAILED.append((name, traceback.format_exc()))
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


# --------------------------------------------------------------------- fixtures
PDB = sorted(x["pdb"] for x in I.targets())[0]
CAND = H.Candidates.from_universe(PDB, k=500)
E = H.score_shipped(CAND)


# --------------------------------------------------------------------- the tests
def t_container_rejects_bad_shapes():
    try:
        H.Candidates(pdb=PDB, n=CAND.n, seq=CAND.seq, fold=CAND.fold,
                     W=np.zeros((10, CAND.n + 1, 3)))
    except ValueError:
        pass
    else:
        raise AssertionError("accepted W with the wrong residue count")
    try:
        H.Candidates(pdb=PDB, n=CAND.n, seq=CAND.seq, fold=CAND.fold,
                     W=np.full((10, CAND.n, 3), np.nan))
    except ValueError:
        pass
    else:
        raise AssertionError("accepted non-finite coordinates")


def t_score_is_native_free():
    """The score must not change when the native label is removed from the container."""
    blind = H.Candidates(pdb=CAND.pdb, n=CAND.n, seq=CAND.seq, fold=CAND.fold,
                         W=CAND.W, PHI=CAND.PHI, PSI=CAND.PSI,
                         nat_ca=None, oracle_rr=None, source="blinded")
    assert np.allclose(H.score_shipped(blind), E, rtol=0, atol=0), \
        "shipped score changed when the native was removed -- it is NOT native-free"
    assert np.allclose(H.score_legacy(blind), H.score_legacy(CAND), rtol=0, atol=0), \
        "legacy score changed when the native was removed"


def t_legacy_needs_torsions():
    co = H.Candidates(pdb=CAND.pdb, n=CAND.n, seq=CAND.seq, fold=CAND.fold, W=CAND.W,
                      nat_ca=CAND.nat_ca, source="coord_only")
    try:
        H.score_legacy(co)
    except ValueError:
        pass
    else:
        raise AssertionError("score_legacy silently produced a number without torsions")


def t_zrank_is_monotone():
    """zrank must not reorder -- the theorem's membership claim rides on the ORDER only."""
    z = H.zrank(E)
    assert np.array_equal(np.argsort(z, kind="stable"), np.argsort(E, kind="stable")), \
        "zrank changed the ordering"


def t_both_arms_see_identical_inputs():
    """THE central invariant. Perturb E; both arms must move together, never one alone."""
    r0 = H.run_both(CAND, E, alpha=0.15, T=0.5, seed=0, n_random_draws=2)
    Ep = E.copy()
    Ep[np.argsort(E)[:200]] += 10.0            # demote the whole classical top-200
    r1 = H.run_both(CAND, Ep, alpha=0.15, T=0.5, seed=0, n_random_draws=2)
    assert r0["c_fixed_rmsd"] != r1["c_fixed_rmsd"], "classical arm ignored the energy change"
    assert r0["q_rmsd"] != r1["q_rmsd"], "quantum arm ignored the energy change"


def t_gate_asserts_subset_reports_equality():
    r = H.run_both(CAND, E, alpha=0.15, T=0.5, seed=0, n_random_draws=2)
    assert r["gate_pass"], "subset-hood -- THE THEOREM -- failed; the harness is mis-wired"
    assert r["gate_equality"], ("equality failed on a trained full-support state; expected the "
                                "empirical regime to hold here")
    assert r["q_pad_in_tail"] == 0, "a padding state reached the readout"


def t_gate_detects_a_mismatched_energy():
    """The gate must FIRE when the arms are deliberately desynchronised."""
    v = np.argsort(E, kind="stable")[[0, 1, 400]]     # a member outside any small prefix
    g = H.gate_set_equality(E, v, m=3)
    assert g["equality"] is False, "gate reported equality for a corrupted set"


def t_matched_control_reproduces_the_quantum_arm():
    """By the theorem, the size-matched classical arm IS the VQE's set -> identical RMSD."""
    r = H.run_both(CAND, E, alpha=0.15, T=0.5, seed=0, n_random_draws=2)
    assert abs(r["q_rmsd"] - r["c_matched_rmsd"]) < 1e-9, (
        f"q={r['q_rmsd']} c_matched={r['c_matched_rmsd']}: the size-matched control must be "
        f"bit-equal to the quantum arm under set equality; a gap means a frame or ordering bug")


def t_readout_is_frame_consistent():
    """The average must not depend on a global rigid motion of the whole pool."""
    rng = np.random.default_rng(0)
    A = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    if np.linalg.det(A) < 0:
        A[:, 0] *= -1
    moved = H.Candidates(pdb=CAND.pdb, n=CAND.n, seq=CAND.seq, fold=CAND.fold,
                         W=CAND.W @ A.T + 7.0, PHI=CAND.PHI, PSI=CAND.PSI,
                         nat_ca=CAND.nat_ca, source="rotated")
    idx = np.arange(30)
    C0, _ = H.readout_uniform(CAND, idx)
    C1, _ = H.readout_uniform(moved, idx)
    assert abs(H._oracle_rmsd(CAND, C0) - H._oracle_rmsd(moved, C1)) < 1e-8, \
        "readout is not invariant to a rigid motion of the pool"


def t_random_null_matches_operator_space():
    """The null must go through the SAME coordinate average, at the SAME m."""
    r = H.arm_random(CAND, m=75, n_draws=8, seed=0)
    assert r["m"] == 75 and r["n_draws"] == 8
    assert r["oracle_rmsd_best"] <= r["oracle_rmsd_mean"] <= r["oracle_rmsd_worst"]
    # and it must be WORSE than the scored selection, or the score carries nothing
    top = H.arm_classical(CAND, E, 75)
    C, _ = H.readout_uniform(CAND, top["cands"])
    assert H._oracle_rmsd(CAND, C) < r["oracle_rmsd_mean"], \
        "the scored top-75 did not beat a random 75-subset on this target"


def t_aggregate_refuses_mixed_bases():
    a = dict(basis="point_cloud", fold=0, pdb="X", q_rmsd=1., c_matched_rmsd=1.,
             c_fixed_rmsd=1., rand_rmsd_mean=1., rand_rmsd_best=1., gate_pass=True,
             q_m=5, q_entropy_bits=1., source="s")
    b = dict(a, basis="built_chain")
    try:
        H.aggregate([a, b])
    except ValueError as e:
        assert "refusing to aggregate mixed bases" in str(e)
    else:
        raise AssertionError("aggregate combined point-cloud with built-chain RMSD")


def t_completion_flag_demands_full_key_set():
    p = H.write("_unittest_incomplete", dict(rows=[{"pdb": "X"}]), required_keys=H.ROW_KEYS)
    with open(p.replace(".json", "") + "_COMPLETE") as fh:
        c = json.load(fh)
    assert c["complete"] is False and c["missing_keys"], \
        "completion flag passed a row missing nearly every key"
    for f in (p, p.replace(".json", "") + "_COMPLETE"):
        os.remove(f)


def t_stats_report_everything_the_brief_demands():
    rng = np.random.default_rng(0)
    a = rng.normal(3.0, .5, 126)
    b = a + rng.normal(0.05, .3, 126)
    folds = rng.integers(0, 5, 126)
    s = H.paired_stats(a, b, folds)
    for k in ("se", "mde", "eff_over_mde", "ci_iid", "ci_fold", "W", "L", "median",
              "worst_degradation", "folds_same_sign", "verdict"):
        assert k in s, f"paired_stats is missing {k}"
    assert abs(s["mde"] - 2.8016 * s["se"]) < 1e-12, "MDE is not 2.8016 x SE"
    assert s["verdict"] in ("NULL", "MEASURED", "TYPE-M ZONE (0.7-1.3x MDE, NOT A RESULT)")


def t_type_m_zone_is_flagged():
    """An effect at ~1.0x its own MDE must NOT be reported as a result."""
    rng = np.random.default_rng(3)
    n = 200
    d = rng.normal(0, 1, n)
    d = d - d.mean()
    se = d.std(ddof=1) / np.sqrt(n)
    a = np.zeros(n) + d + 1.0 * 2.8016 * se        # effect placed at exactly 1.0x MDE
    s = H.paired_stats(a, np.zeros(n))
    assert 0.7 <= s["eff_over_mde"] <= 1.3
    assert "TYPE-M" in s["verdict"], f"1.0x MDE was reported as {s['verdict']}"


def t_merge_tracks_provenance():
    other = H.Candidates.from_uniform_library(PDB, k=200, seed=1)
    u = CAND.merge(other)
    assert u.k == CAND.k + other.k
    Eu = H.score_shipped(u)
    r = H.run_both(u, Eu, alpha=0.15, T=0.5, seed=0, n_random_draws=2)
    assert "c_fixed_share_second_source" in r, "merge provenance not reported"
    assert 0.0 <= r["c_fixed_share_second_source"] <= 1.0


def t_arbitrary_pool_runs_end_to_end():
    """The actual deliverable: a pool that is NOT the shipped retrieval set goes through."""
    gen = H.Candidates.from_uniform_library(PDB, k=333, seed=7, source="fake_generator")
    r = H.run_both(gen, H.score_shipped(gen), alpha=0.15, T=0.5, seed=0, n_random_draws=2)
    assert r["k"] == 333 and r["q_n_qubits"] == 9        # ceil(log2(333)) = 9
    assert r["gate_pass"] and np.isfinite(r["q_rmsd"])
    assert r["source"] == "fake_generator"


def t_gauge_control_available():
    r0 = H.run_both(CAND, E, seed=0, n_random_draws=2, label_perm=False)
    r1 = H.run_both(CAND, E, seed=0, n_random_draws=2, label_perm=True)
    assert r1["q_gauge_permuted"] and not r0["q_gauge_permuted"]
    assert r1["gate_pass"], "gate failed under a permuted label -- the encoding is gauge-broken"


def t_register_ceiling_refuses_rather_than_approximates():
    big = H.Candidates(pdb=CAND.pdb, n=CAND.n, seq=CAND.seq, fold=CAND.fold,
                       W=np.repeat(CAND.W, 40, axis=0)[:20000], nat_ca=CAND.nat_ca)
    try:
        H.arm_vqe(big, np.zeros(big.k))
    except ValueError as e:
        assert "ceiling" in str(e)
    else:
        raise AssertionError("harness silently ran a register above the ceiling")


def main() -> int:
    print(f"harness unit tests  (fixture target {PDB}, k={CAND.k}, n={CAND.n})")
    for nm, fn in sorted(globals().items()):
        if nm.startswith("t_") and callable(fn):
            check(nm, fn)
    print(f"\n{len(FAILED)} failed")
    for nm, tb in FAILED:
        print(f"\n--- {nm} ---\n{tb}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
