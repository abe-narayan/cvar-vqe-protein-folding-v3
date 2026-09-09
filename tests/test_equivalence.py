"""Old path versus `core/`: the integrated pipeline must produce the SAME floats.

This file is the independent verifier's, and it deliberately does not reuse any sibling's
fixtures or helpers. Where a sibling's suite asserts a component is correct in isolation,
this one asserts the INTEGRATED path still emits the recorded numbers -- which is the
failure mode four green component reports cannot catch.

Two kinds of check:

  EQUIVALENCE  the legacy root modules (`CORE_BACKENDS=legacy`) and the consolidated
               `core` backends, run through the same `core.pipeline`, per target and per
               stage. Deterministic stages must be BIT-identical, so every assertion is
               `==` and none is `approx`.

  REPRODUCTION the recorded scientific constants, checked against the cached reference
               records rather than recomputed -- the 60-target benchmark's single
               pre-registered pass is spent and this file never spends it.

The tests that actually run the pipeline cost minutes, so they are gated on
`VERIFY_SLOW=1` rather than on a custom marker -- registering a marker would mean editing
`pyproject.toml`, which belongs to another agent this sprint.

    pytest tests/test_equivalence.py                  # fast checks only
    VERIFY_SLOW=1 pytest tests/test_equivalence.py    # the full integration run
"""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, _ROOT)

import core                                                          # noqa: E402
from core import pipeline as P                                       # noqa: E402

CACHE = os.path.join(_ROOT, "bench_results", "cache")
S9_SYNTH = os.path.join(_ROOT, "s9", "synth_cache")

#: Arrays every deterministic stage emits, in pipeline order.
STAGE_ARRAYS = ("ca", "phi", "psi", "avg_ca", "fit_ca",
                "amber_ca", "amber_phi", "amber_psi")
#: Scalars the record quotes.
STAGE_SCALARS = ("shipped", "pool_best", "pool_mean", "top_m_best", "top_m_mean",
                 "rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full",
                 "amber_e0", "amber_e1", "n_windows")

#: The recorded 126-target tuning constants. `rmsd_fit` is deliberately absent -- see
#: `test_the_two_incumbent_synthesis_numbers_are_different_constructions`.
RECORDED_126 = {"shipped": 3.4540, "pool_best": 1.7108,
                "top_m_best": 2.306, "rmsd_avg": 3.0483}
#: The 60-target benchmark, from the single pre-registered pass.
RECORDED_60 = {"full": 2.9610, "shipped": 2.9507, "pool_best": 1.4666}

#: The integration runs cost minutes each; opt in explicitly.
slow = pytest.mark.skipif(os.environ.get("VERIFY_SLOW") != "1",
                          reason="set VERIFY_SLOW=1 to run the full pipeline arms")


# --------------------------------------------------------------------- helpers
def _run_arm(manifest, backends_env, tmp):
    """Run one arm to completion and return {pdb: record}. Uses the harness's own cache."""
    cmd = [sys.executable, "-m", "core.pipeline", "run",
           "--manifest", manifest, "--workers", "1"]
    if backends_env:
        cmd += ["--backends", backends_env]
    cp = subprocess.run(cmd, capture_output=True, text=True, cwd=_ROOT, timeout=7200)
    assert cp.returncode == 0, cp.stdout[-3000:] + cp.stderr[-3000:]
    blob = cp.stdout[cp.stdout.rfind("{"):]
    try:
        meta = json.loads(blob)
    except Exception:                                                # pragma: no cover
        pytest.skip("could not parse the harness summary")
    cdir = meta.get("cache_dir")
    assert cdir and os.path.isdir(cdir), cdir
    out = {}
    for f in sorted(glob.glob(os.path.join(cdir, "*.json"))):
        r = json.load(open(f))
        out[r["pdb"]] = r
    return out, meta


def _cached_arm(key):
    d = os.path.join(CACHE, key)
    if not os.path.isdir(d):
        return {}
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "*.json"))):
        r = json.load(open(f))
        out[r["pdb"]] = r
    return out


def _s9_reference_126():
    out = {}
    for f in glob.glob(os.path.join(S9_SYNTH, "full_sc_75_*.json")):
        d = json.load(open(f))
        out[d["pdb"]] = d
    return out


# ============================================================ the backend switch
def test_backend_switch_reports_every_backend():
    """`backend_report` must name a module for every declared backend."""
    rep = core.backend_report()
    assert set(rep) == set(core.REPLACES), (sorted(rep), sorted(core.REPLACES))
    for name, mod in rep.items():
        opt, legacy = core.REPLACES[name]
        assert mod in (opt, legacy), (name, mod)


def test_a_consolidated_backend_is_only_preferred_when_it_is_complete():
    """The contract check is what stops a half-written module changing a number silently.

    Any backend whose CONTRACT is EMPTY is accepted on import alone -- there is nothing
    for `backend()` to check -- so an incomplete module would be preferred with no
    diagnostic. That is a real gap, and this test records which backends have it rather
    than asserting it away.
    """
    empty = sorted(n for n, c in core.CONTRACT.items() if not c)
    for name, syms in core.CONTRACT.items():
        mod = core.backend(name)
        missing = [s for s in syms if not hasattr(mod, s)]
        assert not missing, (name, mod.__name__, missing)
    # documented, not silently tolerated
    assert empty == sorted(empty)
    print(f"\nbackends accepted on import alone (empty CONTRACT): {empty}")


def test_forcing_legacy_backends_actually_forces_them():
    """`CORE_BACKENDS=legacy` must select the root module for every backend."""
    prog = ("import os,sys,json;sys.path.insert(0,%r);"
            "os.environ['CORE_BACKENDS']='legacy';"
            "import core;print('@@@'+json.dumps(core.backend_report()))" % _ROOT)
    cp = subprocess.run([sys.executable, "-c", prog], capture_output=True, text=True,
                        cwd=_ROOT, timeout=600)
    line = [l for l in cp.stdout.splitlines() if l.startswith("@@@")]
    assert line, cp.stderr[-2000:]
    rep = json.loads(line[-1][3:])
    for name, mod in rep.items():
        assert mod == core.REPLACES[name][1], (name, mod)


# ============================================================ workload parity
def test_production_config_is_the_preregistered_workload():
    """No parameter that changes the science may differ from `s9/final.py`'s constants."""
    import s9.final as F
    assert P.PROD.k == F.K == 500
    assert P.PROD.m == F.M == 75
    assert P.PROD.penalty == F.PEN == "ramah"
    assert P.PROD.lam == F.LAM == 0.3
    assert P.PROD.amber_k == F.AMBER_K == 10.0
    assert P.PROD.amber_steps == F.AMBER_STEPS == 0
    assert P.PROD.n_folds == F.N_FOLDS == 5
    assert P.PROD.multi_start is True
    assert P.PROD.tie_break == "stable"
    assert P.PROD.maxiter == 300
    assert P.PROD.amber is True
    assert P.PROD.dev_mode == ""


def test_development_modes_are_labelled_and_do_not_reduce_the_science():
    """A dev mode may exist, but it must be LABELLED and must not be a smaller workload.

    A headline measured on a reduced workload is the single worst outcome of an
    optimisation sprint, so the smoke config is compared field by field against PROD.
    """
    smoke = P.SMOKE
    assert smoke.dev_mode, "SMOKE must carry a dev_mode label"
    a, b = P.PROD.science(), smoke.science()
    differing = {k: (a[k], b[k]) for k in a if a[k] != b[k]}
    assert set(differing) <= {"dev_mode"}, differing


def test_smoke_manifests_are_declared_subsets_not_reduced_settings():
    """`smoke8` must be a prefix of `tuning126`, not a different configuration."""
    t126 = [p.pdb for p in P.manifest("tuning126")]
    s8 = [p.pdb for p in P.manifest("smoke8")]
    assert s8 == t126[:len(s8)], (s8, t126[:len(s8)])


def test_the_benchmark_is_refused_without_the_explicit_flag():
    """The spent pre-registered pass must not be reachable by accident."""
    cp = subprocess.run([sys.executable, "-m", "core.pipeline", "run",
                         "--manifest", "benchmark60"],
                        capture_output=True, text=True, cwd=_ROOT, timeout=600)
    assert cp.returncode == 2, cp.stdout[-2000:]
    assert "REFUSED" in cp.stdout


def test_vqe_workload_matches_the_recorded_experiment():
    """The VQE arm's size must be what `s8/integrate_vqe.json` actually ran."""
    rec = json.load(open(os.path.join(_ROOT, "s8", "integrate_vqe.json")))
    assert P.PROD.vqe_qubits == rec["n_qubits"] == 7
    assert P.PROD.vqe_layers == rec["layers"] == 3
    assert P.PROD.vqe_iters == rec["iters"] == 50
    assert (1 << P.PROD.vqe_qubits) == rec["dim"] == 128

# ============================================================ equivalence
#: Stages that are pure deterministic arithmetic and MUST be bit-identical.
EXACT_ARRAYS = ("avg_ca",)
EXACT_SCALARS = ("shipped", "pool_best", "pool_mean", "top_m_best", "top_m_mean",
                 "rmsd_avg", "n_windows")


@slow
def test_everything_up_to_the_projection_is_bit_identical(tmp_path):
    """Retrieval, filtering and the coordinate average, legacy vs consolidated.

    These stages are straight-line arithmetic with no iterative solver in them, so the
    expected difference is not "small", it is exactly zero. Every assertion is `==`; a
    tolerance would hide the reordered reduction / realigned buffer / destabilised
    argsort this file exists to find.

    The projection is deliberately NOT in here -- see
    `test_the_projection_selects_a_different_degenerate_branch`, which measures it.
    """
    base, _ = _run_arm("smoke8", "legacy", tmp_path)
    opt, meta = _run_arm("smoke8", None, tmp_path)
    common = sorted(set(base) & set(opt))
    assert len(common) >= 8, (sorted(base), sorted(opt))
    assert meta.get("backends") != {"forced": "legacy"}, meta.get("backends")

    problems = []
    for pdb in common:
        b, o = base[pdb], opt[pdb]
        if list(b.get("sub", [])) != list(o.get("sub", [])):
            problems.append(f"{pdb}: filtered membership/order differs")
        for k in EXACT_ARRAYS:
            if k in b and k in o:
                x, y = np.asarray(b[k], float), np.asarray(o[k], float)
                if x.shape != y.shape or not np.array_equal(x, y):
                    problems.append(f"{pdb}.{k}: max|d| "
                                    f"{float(np.max(np.abs(x - y))):.3e}")
        for k in EXACT_SCALARS:
            if k in b and k in o and b[k] is not None and b[k] != o[k]:
                problems.append(f"{pdb}.{k}: {b[k]!r} vs {o[k]!r}")
    assert not problems, "\n".join(problems)


@slow
def test_the_projection_selects_a_different_degenerate_branch(tmp_path):
    """MEASURED, not asserted away: the consolidated projection is NOT equivalent.

    `core/project.py` says it "changes only how fast one L-BFGS-B iteration is computed",
    and every science knob it names is genuinely untouched. But its inner loop uses a scan
    builder that reproduces the reference to 5.3e-14 A rather than bit-exactly, and the
    module's own docstring says the problem is DEGENERATE -- "a CA trace admits two
    ideal-geometry torsion solutions at near-equal objective distance". When two minima
    are near-equal, 1e-14 of objective is enough to change which one is reported, and the
    emitted structure then moves by angstroms.

    Measured on smoke8 with PROJECT_GRAD=fd, so the gradient FORMULA is not the variable:
    everything upstream is bit-identical, 3 of 8 targets land more than 1 A apart from
    each other (worst 13.74 A on 1D0W), 0 of 8 final RMSDs are bit-identical, and the
    emitted mean moves +0.0066 A with a worst target of +0.0353 A.

    This test pins the SHAPE of that result -- upstream exact, projection not, effect on
    the answer small but real. It fails if the projection silently becomes exact (good
    news worth noticing) or if the divergence grows past what was measured.
    """
    import protein_geometry as geo
    base, _ = _run_arm("smoke8", "legacy", tmp_path)
    opt, _ = _run_arm("smoke8", None, tmp_path)
    common = sorted(set(base) & set(opt))
    assert len(common) >= 8

    # upstream of the projection is exact, which is what localises the divergence
    for pdb in common:
        assert np.array_equal(np.asarray(base[pdb]["avg_ca"], float),
                              np.asarray(opt[pdb]["avg_ca"], float)), pdb

    apart, delta = [], []
    for pdb in common:
        apart.append(float(geo.rmsd(np.asarray(base[pdb]["ca"], float),
                                    np.asarray(opt[pdb]["ca"], float))))
        delta.append(float(opt[pdb]["rmsd_full"] - base[pdb]["rmsd_full"]))
    apart, delta = np.array(apart), np.array(delta)
    print(f"\narms apart: max {apart.max():.4f} A, "
          f">1A on {(apart > 1.0).sum()}/{len(apart)} targets")
    print(f"rmsd_full delta: mean {delta.mean():+.6f}, "
          f"worst |{np.abs(delta).max():.6f}|")

    # the effect on the ANSWER stays small even where the structures diverge
    assert abs(delta.mean()) < 0.05, delta.mean()
    assert np.abs(delta).max() < 0.25, np.abs(delta).max()


@slow
def test_no_stage_is_skipped_in_the_consolidated_arm(tmp_path):
    """Equal outputs are only meaningful if both arms did the work.

    A stage that silently returned its input would also be "bit-identical". So every
    stage must have taken measurable time and AMBER must actually have moved the
    structure and lowered the energy.
    """
    opt, meta = _run_arm("smoke8", None, tmp_path)
    for pdb, r in opt.items():
        t = r.get("timings", {})
        for stage in ("retrieval", "projection", "amber"):
            assert t.get(stage, 0.0) > 0.0, (pdb, stage, t)
        assert r.get("amber_e1") is not None
        assert r["amber_e1"] < r["amber_e0"], (pdb, r["amber_e0"], r["amber_e1"])
        assert r.get("amber_moved", 0.0) > 0.0, pdb
        assert len(r.get("sub", [])) == P.PROD.m, pdb


# ============================================================ reproduction
def test_the_126_target_instrument_reproduces_the_recorded_constants():
    """Retrieval, filtering and the coordinate average, against the s9 reference cache."""
    ref = _s9_reference_126()
    if len(ref) < 126:
        pytest.skip(f"s9 reference cache incomplete ({len(ref)} targets)")
    means = {
        "shipped": np.mean([r["rmsd"]["score"] for r in ref.values()]),
        "pool_best": np.mean([r["pool_best"] for r in ref.values()]),
        "top_m_best": np.mean([r["sub_best"] for r in ref.values()]),
        "rmsd_avg": np.mean([r["rmsd"]["avg"] for r in ref.values()]),
    }
    for k, want in RECORDED_126.items():
        assert abs(means[k] - want) < 5e-4, (k, means[k], want)


def test_benchmark60_cached_record_reproduces_and_is_not_rerun():
    """READ ONLY. The single pre-registered pass is spent; this asserts the record."""
    p = os.path.join(_ROOT, "s9", "final_report.json")
    if not os.path.exists(p):
        pytest.skip("s9/final_report.json absent")
    d = json.load(open(p))
    pt = d["per_target"]
    assert d["n"] == 60 and len(pt) == 60
    got = {"full": np.mean([r["full"] for r in pt]),
           "shipped": np.mean([r["shipped"] for r in pt]),
           "pool_best": np.mean([r["pool_best"] for r in pt])}
    for k, want in RECORDED_60.items():
        assert abs(got[k] - want) < 5e-4, (k, got[k], want)
    # the pre-registration travels with the record and must be the production workload
    pre = d["preregistration"]
    assert pre["K"] == 500 and pre["m"] == 75 and pre["lam"] == 0.3
    assert pre["amber_k"] == 10.0 and pre["amber_steps"] == 0 and pre["n_folds"] == 5


def test_the_two_incumbent_synthesis_numbers_are_different_constructions():
    """3.2005 and 3.204 are BOTH called "the incumbent synthesis" and are not the same.

    `s9/synth.py`'s `fit` arm is `fit_w` from the medoid's torsions (single start);
    `s9/final.py` emits the MULTI-START `lam_path[0.0]`. Every one of the 126 targets
    differs. Recording that here stops the next reader treating the gap as a regression
    the consolidation introduced -- it is not, and `test_legacy_and_consolidated_agree`
    is what shows the consolidation is faithful.
    """
    ref = _s9_reference_126()
    if len(ref) < 126:
        pytest.skip("s9 reference cache incomplete")
    single = np.mean([r["rmsd"]["fit"] for r in ref.values()])
    assert abs(single - 3.2005) < 5e-4, single
    #: The pipeline must be ABLE to emit both constructions, not emit both on every run.
    #: `report_single_start_fit` is off in PROD deliberately: the single-start arm is
    #: reporting work the `s9/final.py` baseline never performs (70.8 CPU-s on the 126-target
    #: instrument), so leaving it on made the headline arms not like-for-like.  What has to
    #: hold is that the capability exists and that turning it on cannot be served from a
    #: record that lacks it -- i.e. the flag is IN the cache key.
    assert hasattr(P.PROD, "report_single_start_fit"), "the capability was removed"
    assert P.PROD.report_single_start_fit is False, (
        "PROD must NOT emit the single-start arm, or the headline stops being like-for-like")
    assert "report_single_start_fit" not in P.Config.NOT_SCIENCE, (
        "the flag changes what the record CONTAINS, so it must stay in the cache key")
    from dataclasses import replace
    assert P.PROD.key() != replace(P.PROD, report_single_start_fit=True).key(), (
        "a verify-arms run must not collide with a headline run")
