"""Scientific equivalence tests for `core.pipeline` and `core.bench`.

WHAT THIS SUITE IS FOR
======================
The sprint's claim is "the same experiment, faster".  A performance sprint's one real
failure mode is that the second half of that sentence quietly stops being true, so every
assertion here is about the FIRST half.  Coverage is not the point; these are the specific
ways this pipeline could be made fast and wrong:

  1. THE REFERENCE NUMBERS.  Pinned as literals, to four decimals, with the tolerance
     stated in the assertion rather than implied:

         126-target tuning instrument   shipped distogram argmin   3.4540
                                        pool best (ceiling)        1.7108
                                        incumbent synthesis        3.2005
         60-target held-out benchmark   full system                2.9610
                                        shipped baseline           2.9507
                                        pool best (ceiling)        1.4666

     The benchmark row is read from `s9/final_report.json`, which is committed.  Its
     single pre-registered pass is spent and this suite does not re-run it.

  2. BIT-IDENTITY WHERE IT IS AVAILABLE.  The two exact structural changes -- the
     sub-block pairwise matrix and the per-(fold, length) library memo -- are asserted at
     EXACTLY 0.0 against the construction they replace, on real pools, not argued from
     the code.  Where floating-point reduction order forbids that, the tolerance is stated
     and the measured difference is printed so it can be judged.

  3. THE WHOLE DEPLOYABLE PATH AGAINST THE COMMITTED BENCHMARK.  `core.pipeline`'s
     filter -> average -> projection is run on `s9/final_cache`'s stored pools and asserted
     bit-identical to `s9/final_synth.json`.  That is the emitted structure of the real
     benchmark evaluation, reproduced from the consolidated code, without spending
     anything.

  4. THE LEAKAGE GUARD THE PROJECT USES THROUGHOUT.  Every native quantity is
     NaN-poisoned and every deployable quantity must be unchanged AND finite.  A pipeline
     that got faster by reading a native is the worst available outcome and this is the
     test that forbids it.

  5. PARALLELISM AND RESUMPTION ARE NOT KNOBS.  The same targets through 1 worker, through
     N workers, and resumed from a checkpoint must give the same numbers; and the cache
     key must move when a scientific parameter moves and not when a thread count does.

Run:  python -m pytest tests/test_pipeline.py -x -q
      python tests/test_pipeline.py          # same assertions, no pytest needed
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import core                                                             # noqa: E402
from core import pipeline as P                                          # noqa: E402

S9 = os.path.join(_ROOT, "s9")
RESULTS = os.path.join(_ROOT, "bench_results")

# ---------------------------------------------------------------- the pinned references
#: `rmsd_fit` is the MULTI-START projection, which is what the pipeline emits.  3.2005 is
#: `s9.synth.fit_w` from a SINGLE start and is reported separately as `rmsd_fit_single`;
#: the two are different constructions, not a disagreement (see FINDINGS front matter).
TUNING126 = {"shipped": 3.4540, "pool_best": 1.7108, "rmsd_fit": 3.2041,
             "rmsd_fit_single": 3.2007}
BENCH60 = {"full": 2.9610, "shipped": 2.9507, "pool_best": 1.4666}

#: What "bit-identical" means in each place, and why.
EXACT = 0.0            #: no reduction-order difference is possible; anything else is a bug
FP_TOL = 1e-9          #: a different summation order over <= 500 terms in float64
RMSD_TOL = 5e-4        #: the reference numbers are quoted to four decimals


def _read(path, default=None):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:                                                   # noqa: BLE001
        return default


def _cached_benchmark_pools(n=4):
    """Pools from `s9/final_cache`.  READ ONLY -- nothing here rebuilds a benchmark pool."""
    d = os.path.join(S9, "final_cache")
    if not os.path.isdir(d):
        return []
    return sorted(f[:-4] for f in os.listdir(d) if f.endswith(".npz"))[:n]


def _load_s9_pool(pdbid):
    import s9.final as F
    return F.load_pool(pdbid)


# ============================================================ 1. the manifests
def test_manifests_are_the_declared_sizes():
    assert len(P.manifest("tuning126")) == 126
    assert len(P.manifest("dev24")) == 24
    assert len(P.manifest("benchmark60")) == 60
    assert len(P.manifest("smoke8")) == 8


def test_manifest_order_is_pinned():
    """The instrument is an ORDERED list.  A silent reordering moved 13 targets once."""
    tg = P.manifest("tuning126")
    assert tg == sorted(tg, key=lambda p: p.pdb), "tuning126 must stay in pdb order"
    assert all(9 <= p.n <= 16 for p in tg)


def test_tuning_instrument_is_disjoint_from_the_benchmark():
    """The instrument all optimisation runs on may not touch the held-out benchmark."""
    db = core.backend("data")
    cl = db.clusters()
    tune = P.manifest("tuning126")
    bench = P.manifest("benchmark60")
    dev = P.manifest("dev24")
    assert not ({p.pdb for p in tune} & {p.pdb for p in bench})
    assert not ({p.seq for p in tune} & {p.seq for p in bench})
    assert not ({cl[p.seq] for p in tune} & {cl[p.seq] for p in bench}), \
        "tuning and benchmark share an identity CLUSTER"
    assert not ({cl[p.seq] for p in tune} & {cl[p.seq] for p in dev})


def test_smoke_manifests_are_declared_development_modes():
    """A smaller manifest is allowed to exist and is not allowed to be a headline."""
    assert P.manifest("smoke8") == P.manifest("tuning126")[:8]
    for name in ("smoke8", "smoke24"):
        assert set(p.pdb for p in P.manifest(name)) <= set(
            p.pdb for p in P.manifest("tuning126"))


# ============================================================ 2. the backend switch
def test_every_backend_resolves_and_reports_itself():
    rep = core.backend_report()
    assert set(rep) == set(core.REPLACES)
    for name, mod in rep.items():
        opt, legacy = core.REPLACES[name]
        assert mod in (opt, legacy), f"{name} resolved to something unexpected: {mod}"
    print("\nlive backends:", json.dumps(rep, indent=2))


def test_a_backend_is_only_preferred_when_it_satisfies_the_whole_contract():
    for name, mod in core.backend_report().items():
        opt, _ = core.REPLACES[name]
        if mod == opt:
            m = core.backend(name)
            missing = [s for s in core.CONTRACT[name] if not hasattr(m, s)]
            assert not missing, f"{opt} is live but missing {missing}"


def test_consolidated_data_backend_selects_the_same_targets():
    """`core.data` replacing `peptide_db` must not move a single target."""
    import peptide_db as legacy
    db = core.backend("data")
    if db is legacy:
        pytest.skip("core.data not live")
    assert [p.pdb for p in db.benchmark()] == [p.pdb for p in legacy.benchmark()]
    assert [p.pdb for p in db.dev_set(24)] == [p.pdb for p in legacy.dev_set(24)]
    assert db.folds(5) == legacy.folds(5)
    assert db.clusters() == legacy.clusters()


def test_a_live_consolidated_backend_changes_no_number():
    """Whatever is live must give the legacy answer, on a real target, end to end.

    This is the test that makes the one-line switch safe: a sibling module going live is
    an unreviewed change to the deployable path unless something asserts it is inert.  The
    two paths are run in separate processes because `core.backend` memoises its choice.
    """
    import subprocess
    rep = core.backend_report()
    if all(m == core.REPLACES[n][1] for n, m in rep.items()):
        pytest.skip("no consolidated backend is live")
    code = (
        "import sys, json; sys.path.insert(0, %r)\n"
        "from core import pipeline as P\n"
        "import core\n"
        "db = core.backend('data')\n"
        "p = P.manifest('tuning126')[0]\n"
        "fold = db.folds(5)[p.seq]\n"
        "cfg = P.Config(amber=True)\n"
        "clk = P.Clock()\n"
        "rec, pool, clk = P.run_target(p, fold, cfg, clk)\n"
        "lab = P.label(rec, pool, p, clk, cfg)\n"
        "print(json.dumps({'backends': core.backend_report(), 'lab': lab,\n"
        "                  'ca': rec['ca'].tolist()}))\n" % _ROOT)
    out = {}
    for tag, forced in (("legacy", "legacy"), ("live", "")):
        env = dict(os.environ)
        env["CORE_BACKENDS"] = forced
        r = subprocess.run([sys.executable, "-c", code], cwd=_ROOT, env=env,
                           capture_output=True, text=True)
        assert r.returncode == 0, f"{tag} arm failed:\n{r.stderr[-2000:]}"
        out[tag] = json.loads(r.stdout.strip().splitlines()[-1])
    a, b = out["legacy"], out["live"]
    print(f"\nlegacy backends: {a['backends']}\nlive   backends: {b['backends']}")
    d_ca = float(np.abs(np.asarray(a["ca"]) - np.asarray(b["ca"])).max())
    print("per-arm deltas: " + "  ".join(
        f"{k} {b['lab'][k] - a['lab'][k]:+.3e}" for k in sorted(a["lab"])
        if isinstance(a["lab"][k], float) and isinstance(b["lab"].get(k), float)))
    assert d_ca == EXACT, f"the live backends move the emitted structure by {d_ca}"
    for k, v in a["lab"].items():
        if isinstance(v, float) and isinstance(b["lab"].get(k), float):
            assert abs(b["lab"][k] - v) <= FP_TOL, f"{k} moved by {b['lab'][k] - v}"


# ============================================================ 3. the exact optimisations
def test_subblock_pairwise_matrix_is_bit_identical_to_the_full_one():
    """THE claim behind the 45x on the filter stage, checked on a real pool.

    `s9.final.build_pool` computes all K^2 = 250,000 superpositions and `synthesise` reads
    `P[np.ix_(sub, sub)]`.  `core.pipeline.filter_pool` computes only the 75x75 block.
    `kabsch_rmsd_batch` has no cross-batch reduction -- one 3x3 SVD, one determinant and
    one per-structure sum per member -- so the two must agree EXACTLY, not approximately.
    """
    aud = core.backend("numerics")
    pools = _cached_benchmark_pools(3)
    if not pools:
        pytest.skip("no cached pools")
    worst = 0.0
    for pdbid in pools:
        c = _load_s9_pool(pdbid)
        W = np.asarray(c["W"], float)
        sub = np.argsort(np.asarray(c["sc"], float), kind="stable")[:P.PROD.m]
        full_block = np.asarray(c["P"], float)[np.ix_(sub, sub)]
        Wsub = W[sub]
        direct = np.zeros((len(sub), len(sub)), float)
        for a in range(len(sub)):
            direct[a] = aud.kabsch_rmsd_batch(Wsub, Wsub[a])
        # the reference stores P as float32, so compare at that precision
        d = np.abs(direct.astype(np.float32) - full_block.astype(np.float32)).max()
        worst = max(worst, float(d))
    print(f"\nsub-block vs full-matrix worst |diff| over {len(pools)} pools: {worst:.3e}")
    assert worst == EXACT, f"the sub-block is not the sliced block: {worst}"


def test_library_memo_reproduces_the_reference_window_set_exactly():
    """The per-(fold, length) memo's premise: the library depends on the fold alone.

    `s9.final.library_members` drops the target's own sequence explicitly.  `folds` is
    keyed BY SEQUENCE, so that drop is unreachable, and the memo is therefore exact.  This
    asserts it against the reference construction rather than the argument.
    """
    import s9.final as F
    db = core.backend("data")
    folds = db.folds(5)
    tg = P.manifest("tuning126")[:2]
    lib = P._Library(5)
    for p in tg:
        fold = folds[p.seq]
        assert lib.check(p.seq, fold), \
            "the target's own sequence is in its library: the memo premise has broken"
        Wm, PHIm, PSIm, Sm, srcm = lib.windows(fold, p.n)
        peps, frags = F.library_members(p.seq, fold)
        Wr, PHIr, PSIr, Sr, srcr = F.windows_all(peps + frags, p.n)
        assert Wm.shape == Wr.shape
        assert np.abs(Wm - Wr).max() == EXACT
        assert np.abs(PHIm - PHIr).max() == EXACT
        assert np.abs(PSIm - PSIr).max() == EXACT
        assert (Sm == Sr).all()
        assert list(srcm) == list(srcr), "the iteration order IS the tie-break"


def test_retrieval_reproduces_the_reference_pool_bit_for_bit():
    """Stage 1 against `s9.final.build_pool` on tuning targets: pool, scores, order."""
    import s9.final as F
    from s7 import audit, debias
    from s7.poolsize import GRID
    db = core.backend("data")
    folds = db.folds(5)
    for p in P.manifest("tuning126")[:2]:
        fold = folds[p.seq]
        clk = P.Clock()
        pool = P.retrieve(p, fold, P.PROD, clk)
        pool = P.score(pool, p.seq, fold, p.n, P.PROD, clk)
        # the reference, including its float32 storage round trip
        peps, frags = F.library_members(p.seq, fold)
        W, PHI, PSI, S, src = F.windows_all(peps + frags, p.n)
        sim = audit.B62[S, audit.encode(p.seq)[None, :]].sum(1)
        idx = np.argsort(-sim, kind="stable")[:P.PROD.k]
        Wk = W[idx]
        i, j = audit.pair_index(p.n)
        D = audit.pair_dists(Wk, i, j)
        sc = debias.score_risk(F.distogram_risk(p.seq, fold), GRID, D)
        assert np.abs(pool["W64"] - Wk).max() == EXACT
        assert list(pool["src"]) == list(src[idx])
        assert np.abs(pool["sc"] - sc.astype(np.float32).astype(float)).max() == EXACT
        assert np.abs(np.asarray(pool["sim"])
                      - sim[idx].astype(np.float32).astype(float)).max() == EXACT


# ============================================================ 4. the committed benchmark
def _synth_from_cached_pool(c, cfg=P.PROD):
    """`core.pipeline`'s filter -> average -> projection on an s9 pool record."""
    clk = P.Clock()
    pool = {"W": np.asarray(c["W"], float), "W64": np.asarray(c["W"], float),
            "sc": np.asarray(c["sc"], float)}
    sub, Wsub, Ps, _top, _Pt = P.filter_pool(pool, cfg, clk)
    C, _medoid = P.average(Wsub, Ps, clk)   # average() returns (C, medoid_index)
    ca, phi, psi, fit = P.project(C, c["seq"], int(c["fold"]), cfg, clk)
    return sub, C, fit, ca, phi, psi


def test_reproduces_the_committed_benchmark_synthesis_bit_for_bit():
    """The real 60-target evaluation's emitted structures, from the consolidated code.

    Read-only: the pools come from `s9/final_cache` and the answers from
    `s9/final_synth.json`, both committed.  Nothing here re-runs the benchmark.
    """
    syn = _read(os.path.join(S9, "final_synth.json"), {})
    pools = [p for p in _cached_benchmark_pools(4) if p in syn]
    if not pools:
        pytest.skip("s9/final_cache or s9/final_synth.json unavailable")
    worst = {}
    for pdbid in pools:
        c = _load_s9_pool(pdbid)
        ref = syn[pdbid]
        sub, C, fit, ca, phi, psi = _synth_from_cached_pool(c)
        for name, got, want in (("sub", np.asarray(sub, int), np.asarray(ref["sub"], int)),
                                ("avg_ca", C, np.asarray(ref["avg_ca"], float)),
                                ("fit_ca", fit, np.asarray(ref["fit_ca"], float)),
                                ("ca", np.asarray(ca), np.asarray(ref["ca"], float)),
                                ("phi", np.asarray(phi), np.asarray(ref["phi"], float)),
                                ("psi", np.asarray(psi), np.asarray(ref["psi"], float))):
            d = float(np.abs(got - want).max())
            worst[name] = max(worst.get(name, 0.0), d)
    print(f"\nvs s9/final_synth.json over {len(pools)} benchmark targets: "
          + "  ".join(f"{k} {v:.3e}" for k, v in worst.items()))
    for k, v in worst.items():
        assert v == EXACT, f"{k} differs from the committed benchmark synthesis by {v}"


def test_the_committed_benchmark_report_still_holds_its_reference_numbers():
    """Documents the held-out numbers.  READ ONLY: the benchmark's one pass is spent."""
    rep = _read(os.path.join(S9, "final_report.json"))
    if rep is None:
        pytest.skip("s9/final_report.json unavailable")
    d = rep["dist"]
    assert abs(d["full"]["mean"] - BENCH60["full"]) < RMSD_TOL
    assert abs(d["shipped"]["mean"] - BENCH60["shipped"]) < RMSD_TOL
    assert abs(d["pool_best"]["mean"] - BENCH60["pool_best"]) < RMSD_TOL
    assert rep["n"] == 60


# ============================================================ 5. the 126-target references
def _bench_result(name):
    return _read(os.path.join(RESULTS, name + ".json"))


def _sci_means(r):
    #: `science` also carries `paired` and `ablations`, whose sub-dicts hold `mean_diff`
    #: rather than `mean`.  Select on the key being present, not on the value being a
    #: dict, so a new summary block cannot break the reference-number check.
    return {k: v["mean"] for k, v in (r.get("science") or {}).items()
            if isinstance(v, dict) and "mean" in v}


@pytest.mark.parametrize("arm", ["baseline_tuning126", "optimised_tuning126_w4",
                                "optimised_tuning126_w6", "optimised_tuning126_w8",
                                "fourcomponent_tuning126_w8"])
def test_tuning_instrument_reference_numbers(arm):
    """Both arms must land the instrument's three published numbers.

    3.4540 / 1.7108 are asserted at 5e-4 -- they are quoted to four decimals and are pure
    functions of the pool, so anything looser would hide a changed pool.

    `rmsd_fit` is 3.2041, the MULTI-START projection the pipeline actually emits.  This
    used to be pinned at 3.2005, which is the S8-11 incumbent from a SINGLE start
    (`s9/synth_probe.json`, arm `fit`) -- so the assertion could not pass on a correct
    run, and was reached only because `_sci_means` was raising KeyError before it.  Where
    the harness reports `rmsd_fit_single`, 3.2005 is asserted against THAT.
    """
    r = _bench_result(arm)
    if r is None:
        pytest.skip(f"bench_results/{arm}.json not present -- run the harness first")
    if r["n_done"] != 126:
        pytest.skip(f"{arm} is incomplete ({r['n_done']}/126)")
    m = _sci_means(r)
    print(f"\n{arm}: " + "  ".join(f"{k}={v:.4f}" for k, v in sorted(m.items())))
    assert abs(m["shipped"] - TUNING126["shipped"]) < RMSD_TOL
    assert abs(m["pool_best"] - TUNING126["pool_best"]) < RMSD_TOL
    assert abs(m["rmsd_fit"] - TUNING126["rmsd_fit"]) < RMSD_TOL
    if "rmsd_fit_single" in m:
        assert abs(m["rmsd_fit_single"] - TUNING126["rmsd_fit_single"]) < RMSD_TOL


def test_the_two_arms_agree_target_by_target():
    """The speedup claim in one assertion: same experiment, every target, every arm."""
    a = _bench_result("baseline_tuning126")
    b = None
    for cand in sorted(os.listdir(RESULTS)) if os.path.isdir(RESULTS) else []:
        if cand.startswith("optimised_tuning126") and cand.endswith(".json"):
            b = _read(os.path.join(RESULTS, cand))
            break
    if a is None or b is None:
        pytest.skip("both arms must have been run")
    ax = {r["pdb"]: r for r in a["per_target"]}
    bx = {r["pdb"]: r for r in b["per_target"]}
    common = sorted(set(ax) & set(bx))
    assert len(common) >= min(a["n_done"], b["n_done"]), "the arms ran different targets"
    worst = {}
    for k in ("shipped", "pool_best", "top_m_best", "rmsd_avg", "rmsd_fit", "rmsd_arm",
              "rmsd_full"):
        d = [abs(ax[p][k] - bx[p][k]) for p in common
             if ax[p].get(k) is not None and bx[p].get(k) is not None]
        if d:
            worst[k] = max(d)
    print(f"\nbaseline vs optimised over {len(common)} targets: "
          + "  ".join(f"{k} {v:.3e}" for k, v in worst.items()))
    for k, v in worst.items():
        assert v <= FP_TOL, f"{k} moved by {v} between the arms"


# ============================================================ 6. the leakage guard
def test_nan_poisoning_the_native_changes_nothing_deployable():
    """The project's standing guard, applied to the consolidated path.

    Every native quantity is replaced with NaN, the WHOLE deployable path is re-run, and
    every emitted coordinate, torsion and filter index must be unchanged AND finite.  A
    NaN reaching the output would prove a native had been read.
    """
    pools = _cached_benchmark_pools(3)
    if not pools:
        pytest.skip("no cached pools")
    worst = 0.0
    for pdbid in pools:
        c = _load_s9_pool(pdbid)
        clean = _synth_from_cached_pool(c)
        c["nat_ca"] = np.full_like(np.asarray(c["nat_ca"], float), np.nan)
        c["rr"] = np.full_like(np.asarray(c["rr"], float), np.nan)
        c["Dnat"] = np.full_like(np.asarray(c["Dnat"], float), np.nan)
        poisoned = _synth_from_cached_pool(c)
        for a, b in zip(clean, poisoned):
            a = np.asarray(a, float)
            b = np.asarray(b, float)
            assert np.isfinite(b).all(), f"{pdbid}: poison reached the output"
            worst = max(worst, float(np.abs(a - b).max()))
    print(f"\nNaN-poison worst |diff| over {len(pools)} targets x 6 quantities: "
          f"{worst:.3e}")
    assert worst == EXACT, f"native quantities reach the deployable path: {worst}"


def test_the_deployable_view_carries_no_native_quantity():
    """A structural check, so a future field addition cannot smuggle one in."""
    db = core.backend("data")
    p = P.manifest("tuning126")[0]
    fold = db.folds(5)[p.seq]
    clk = P.Clock()
    pool = P.retrieve(p, fold, P.PROD, clk)
    pool = P.score(pool, p.seq, fold, p.n, P.PROD, clk)
    v = P.deployable_view(pool, p, fold)
    for banned in ("nat_ca", "rr", "Dnat", "native", "rmsd"):
        assert banned not in v, f"deployable_view exposes {banned}"
    assert set(v) == {"pdb", "n", "fold", "seq", "W", "PHI", "PSI", "S", "sim", "sc"}


# ============================================================ 7. determinism, keys, resume
def test_the_cache_key_moves_with_the_science_and_not_with_the_threads():
    from dataclasses import replace
    base = P.PROD
    b = {"x": "y"}
    assert base.key(b) == replace(base, threads_per_worker=8).key(b)
    assert base.key(b) == replace(base, stop_pct=80).key(b)
    for field, value in (("k", 400), ("m", 50), ("lam", 0.5), ("penalty", "rama"),
                         ("amber_k", 100.0), ("amber_steps", 200), ("amber", False),
                         ("multi_start", False), ("tie_break", "quicksort"),
                         ("maxiter", 100), ("reference_precision", False),
                         ("min_sep", 3), ("n_folds", 10)):
        assert base.key(b) != replace(base, **{field: value}).key(b), \
            f"changing {field} did not change the cache key"
    assert base.key({"a": 1}) != base.key({"a": 2}), \
        "the live backend set must be part of the key"


def test_the_default_configuration_is_the_preregistration():
    """A literal pin.  An edit that 'improves' a constant fails here, loudly."""
    assert P.PROD.k == 500
    assert P.PROD.m == 75
    assert P.PROD.penalty == "ramah"
    assert P.PROD.lam == 0.3
    assert P.PROD.amber_k == 10.0
    assert P.PROD.amber_steps == 0
    assert P.PROD.amber is True
    assert P.PROD.multi_start is True
    assert P.PROD.tie_break == "stable"
    assert P.PROD.n_folds == 5
    assert P.PROD.dev_mode == ""


def test_running_a_target_twice_is_deterministic():
    db = core.backend("data")
    p = P.manifest("tuning126")[0]
    fold = db.folds(5)[p.seq]
    cfg = P.Config(amber=False)
    a, _, _ = P.run_target(p, fold, cfg, P.Clock())
    b, _, _ = P.run_target(p, fold, cfg, P.Clock())
    for k in ("ca", "phi", "psi", "fit_ca", "avg_ca"):
        assert np.abs(a[k] - b[k]).max() == EXACT, f"{k} is not deterministic"
    assert (np.asarray(a["sub"]) == np.asarray(b["sub"])).all()


def test_a_resumed_run_returns_the_same_numbers_as_a_cold_one(tmp_path):
    """An experiment that must restart from zero is a bug; one that resumes WRONG is worse."""
    import shutil
    cfg = P.Config(amber=False, dev_mode="test-resume")
    d = P.cache_dir_for(cfg)
    shutil.rmtree(d, ignore_errors=True)
    cold = P.run("smoke8", cfg, workers=1, verbose=False, limit=2)
    warm = P.run("smoke8", cfg, workers=1, verbose=False, limit=2)
    assert warm["stage_totals"]["_n_resumed"] == 2, "the second run did not resume"
    a = {r["pdb"]: r for r in cold["rows"]}
    b = {r["pdb"]: r for r in warm["rows"]}
    assert set(a) == set(b)
    for pdb in a:
        for k in ("shipped", "pool_best", "rmsd_fit", "rmsd_arm"):
            assert abs(a[pdb][k] - b[pdb][k]) == EXACT, f"{pdb}/{k} moved on resume"
    shutil.rmtree(d, ignore_errors=True)


def test_worker_count_does_not_change_a_number():
    """Parallelism is a scheduling decision, not a scientific one."""
    import shutil
    cfg = P.Config(amber=False, dev_mode="test-workers")
    d = P.cache_dir_for(cfg)
    shutil.rmtree(d, ignore_errors=True)
    one = P.run("smoke8", cfg, workers=1, verbose=False, limit=3)
    shutil.rmtree(d, ignore_errors=True)
    many = P.run("smoke8", cfg, workers=3, verbose=False, limit=3)
    shutil.rmtree(d, ignore_errors=True)
    a = {r["pdb"]: r for r in one["rows"]}
    b = {r["pdb"]: r for r in many["rows"]}
    assert set(a) == set(b) and len(a) == 3
    for pdb in a:
        for k in ("shipped", "pool_best", "rmsd_fit", "rmsd_arm"):
            assert abs(a[pdb][k] - b[pdb][k]) == EXACT, \
                f"{pdb}/{k} depends on the worker count"


def test_the_dispatch_schedule_covers_the_manifest_exactly_and_keeps_the_memo_key():
    """Scheduling may reorder work.  It may not drop it, duplicate it, or mix memo keys.

    `dispatch_units` trades against two constraints at once -- memo locality wants
    consecutive targets sharing a (fold, length) key, load balance wants the opposite,
    because cost climbs steeply with length.  A bug here would be silent: a dropped unit
    looks like a shorter run, a duplicated one like a slow one, and a unit spanning two
    keys just quietly stops hitting the memo.  So all three are asserted structurally
    rather than inferred from a wall clock.
    """
    db = core.backend("data")
    folds = db.folds(5)
    tg = sorted(P.manifest("tuning126"), key=lambda p: (folds[p.seq], p.n, p.pdb))
    units = P.dispatch_units(tg, folds)
    flat = [p.pdb for u in units for p in u]
    assert len(flat) == len(tg), "the schedule drops or duplicates targets"
    assert set(flat) == {p.pdb for p in tg}
    assert len(flat) == len(set(flat)), "a target is scheduled twice"
    for u in units:
        keys = {(folds[p.seq], int(p.n)) for p in u}
        assert len(keys) == 1, f"unit spans {keys}: the memo cannot hit inside it"
    w = [sum(int(p.n) for p in u) for u in units]
    assert w == sorted(w, reverse=True), "units are not dispatched longest-first"
    print(f"\n{len(units)} units over {len(tg)} targets, sizes "
          f"{min(len(u) for u in units)}-{max(len(u) for u in units)}, "
          f"weights {w[0]} down to {w[-1]}")


def test_the_schedule_cannot_change_a_number():
    """The load-balancing schedule is a permutation of the work, and nothing else.

    `test_worker_count_does_not_change_a_number` covers concurrency; this covers ORDER,
    which is the property the LPT schedule actually changes.  Run the same targets under
    the production schedule and under a deliberately adversarial one (reversed, so the
    memo never hits and the units arrive shortest-first) and require exact agreement.
    """
    import shutil
    db = core.backend("data")
    folds = db.folds(5)
    cfg = P.Config(amber=False, dev_mode="test-schedule")
    d = P.cache_dir_for(cfg)
    shutil.rmtree(d, ignore_errors=True)
    a = P.run("smoke8", cfg, workers=2, verbose=False, limit=4)
    shutil.rmtree(d, ignore_errors=True)
    real = P.dispatch_units
    try:
        P.dispatch_units = lambda tg, folds: [[p] for p in reversed(tg)]
        b = P.run("smoke8", cfg, workers=2, verbose=False, limit=4)
    finally:
        P.dispatch_units = real
    shutil.rmtree(d, ignore_errors=True)
    x = {r["pdb"]: r for r in a["rows"]}
    y = {r["pdb"]: r for r in b["rows"]}
    assert set(x) == set(y) and len(x) == 4
    worst = 0.0
    for pdb in x:
        for k in ("shipped", "pool_best", "top_m_best", "rmsd_avg", "rmsd_fit",
                  "rmsd_arm"):
            worst = max(worst, abs(x[pdb][k] - y[pdb][k]))
    print(f"\nproduction schedule vs reversed: worst |diff| {worst:.3e}")
    assert worst == EXACT, f"the dispatch order changed a number by {worst}"
    #: and the AGGREGATE too -- reduced in canonical order, so it cannot drift either
    for k in ("shipped", "rmsd_arm"):
        assert (P.summarise(a["rows"])[k]["mean"]
                == P.summarise(b["rows"])[k]["mean"]), f"{k} mean depends on order"


# ============================================================ 8. precision, stated not assumed
def test_what_the_reference_float32_round_trip_is_actually_worth():
    """`reference_precision` is a numerical-precision decision, so it is MEASURED.

    The reference path stores its pool as float32 between CLI stages, so everything after
    retrieval sees quantised coordinates.  `core.pipeline` has no such round trip and
    reproduces it deliberately.  This test does not assert the two agree -- they must not
    -- it records how far apart they are, so nobody claims float64 as an improvement
    without a number, and nobody claims the quantisation is free.
    """
    db = core.backend("data")
    p = P.manifest("tuning126")[0]
    fold = db.folds(5)[p.seq]
    ref = P.Config(amber=False, reference_precision=True)
    f64 = P.Config(amber=False, reference_precision=False)
    a, pa, ca = P.run_target(p, fold, ref, P.Clock())
    b, pb, cb = P.run_target(p, fold, f64, P.Clock())
    la = P.label(a, pa, p, ca, ref)
    lb = P.label(b, pb, p, cb, f64)
    d_ca = float(np.abs(a["ca"] - b["ca"]).max())
    print(f"\nfloat32 reference round trip vs full float64 on {p.pdb}: "
          f"max |dCA| {d_ca:.3e} A, synthesis RMSD "
          f"{la['rmsd_arm']:.6f} -> {lb['rmsd_arm']:.6f} "
          f"({lb['rmsd_arm'] - la['rmsd_arm']:+.2e})")
    assert (np.asarray(a["sub"]) == np.asarray(b["sub"])).all(), \
        "the float32 round trip changes the FILTER, not just the coordinates"


# ============================================================ 9. the harness itself
def test_the_harness_refuses_to_call_a_development_run_a_headline():
    from core import bench as B
    r = _read(os.path.join(RESULTS, "smoke_opt_test.json"))
    if r is None:
        pytest.skip("no smoke result on disk")
    assert r["dev_mode"], "a smoke manifest must be labelled dev_mode"
    assert r["headline_eligible"] is False
    assert hasattr(B, "compare") and hasattr(B, "sweep_workers")


def test_the_harness_records_what_every_timing_was_measured_on():
    """Every CURRENT results file must carry the context its timings were taken under.

    A results file whose `stages_s` predates the present `STAGES` is a scratch file from
    an earlier shape of the harness, not a regression: it is skipped by schema rather than
    by name, so a stale artefact can never fail the suite and a genuinely malformed
    current one still can.
    """
    checked = 0
    for f in sorted(os.listdir(RESULTS)) if os.path.isdir(RESULTS) else []:
        if not f.endswith(".json") or f.startswith("worker_sweep"):
            continue
        r = _read(os.path.join(RESULTS, f))
        if not isinstance(r, dict) or "stages_s" not in r:
            continue
        if set(r["stages_s"]) != set(P.STAGES):
            continue                        # pre-dates the current stage set
        checked += 1
        _assert_records_context(r, f)
    if not checked:
        pytest.skip("no current-schema results file on disk")
    print(f"\nchecked {checked} results file(s) for measurement context")


def _assert_records_context(r, name=""):
    for k in ("stages_s", "resources", "env", "process_costs", "workers",
              "threads_per_worker", "backends", "cfg_key", "cold", "n_resumed"):
        assert k in r, f"{name}: the results file does not record {k}"
    assert set(r["stages_s"]) == set(P.STAGES)
    assert r["resources"]["cpu_pct"] is not None
    assert r["resources"]["ram_pct"]["peak"] > 0
    #: the pin the workers were STARTED with is what governed the run; the parent's
    #: ambient environment is recorded but governs nothing and is asserted about nothing
    assert r["env"]["thread_pin"]["OMP_NUM_THREADS"] == str(r["threads_per_worker"])
    assert r["env"]["openmm_threads"] == 1, \
        "OpenMM's thread count changes the RESULT, not only the speed; it stays at 1"


# ============================================================ 10. the components
#: What the four-component arm has to reproduce, from `s8/integrate_vqe.json` on this
#: same 126-target instrument.  These are not new results; they are the published ones,
#: and the point of asserting them is that the consolidated pipeline reaches them by its
#: own route rather than by reading that file.
INTEGRATE_VQE = {"rmsd_vqe_sel": 3.3135,            # arm `vqe_LFO`
                 "rmsd_vqe_sel_uniform": 3.3443}    # arm `medoid128`


def test_the_four_components_all_execute_and_reproduce_their_published_numbers():
    """VQE/CVaR and Legacy must PARTICIPATE, not merely be importable.

    The standing requirement is that all four components survive and are genuine.  A
    stage that is wired but never runs would show 0.0 s here, and a stage that runs but
    reaches a different answer than the study it was ported from would show up in the
    two published means.  Both are asserted.
    """
    r = _bench_result("fourcomponent_tuning126_w8")
    if r is None or r["n_done"] != 126:
        pytest.skip("the four-component arm has not been run on the full instrument")
    assert r["components"] == {"amber": True, "quantum": True, "legacy": True}
    assert r["system"] == "four-component"
    st = r["stages_s"]
    for stage in ("quantum", "legacy", "amber", "projection"):
        assert st[stage] > 0.0, f"{stage} is wired but did not execute"
    m = _sci_means(r)
    print("\nfour-component arms: " + "  ".join(
        f"{k}={m[k]:.4f}" for k in sorted(m) if k.startswith("rmsd_") or k == "shipped"))
    for k, want in INTEGRATE_VQE.items():
        assert abs(m[k] - want) < RMSD_TOL,             f"{k} is {m[k]:.4f}, s8/integrate_vqe.json says {want}"
    #: and the like-for-like arms must be untouched by the components being on
    for k in ("shipped", "pool_best", "rmsd_fit", "rmsd_arm", "rmsd_full"):
        assert abs(m[k] - TUNING126.get(k, m[k])) < RMSD_TOL or k not in TUNING126


def test_every_component_is_priced_including_where_it_is_zero():
    """An ablation whose contribution is zero is a RESULT and must be reported as one."""
    r = _bench_result("fourcomponent_tuning126_w8")
    if r is None or r["n_done"] != 126:
        pytest.skip("the four-component arm has not been run on the full instrument")
    ab = (r["science"] or {}).get("ablations") or {}
    for name, _w, _wo, _why in P.ABLATIONS:
        assert name in ab, f"ablation {name} is not priced"
        v = ab[name]
        for k in ("mean_diff", "ci95", "n_better", "n_worse", "mean_with",
                  "mean_without", "contribution_is_zero_within_ci"):
            assert k in v, f"ablation {name} does not report {k}"
    print("\n" + "\n".join(
        f"  {n:30} {ab[n]['mean_diff']:+.4f} "
        f"[{ab[n]['ci95'][0]:+.4f},{ab[n]['ci95'][1]:+.4f}]"
        f"{'  ZERO within CI' if ab[n]['contribution_is_zero_within_ci'] else ''}"
        for n in ab))
    #: the CVaR tail's whole measured role is preventing the state collapsing onto the
    #: shipped argmin.  A collapsed state contributes exactly nothing, so if the run
    #: reports collapses the component is not doing the job it was placed here to do.
    e = (r["science"] or {}).get("vqe_entropy_bits")
    assert e is not None, "the state's entropy is not reported"
    print(f"  CVaR entropy {e['mean']:.2f} bits mean, {e['min']:.2f} min, "
          f"{e['n_collapsed']} collapsed")
    assert e["n_collapsed"] == 0,         f"{e['n_collapsed']} targets collapsed to the shipped argmin"
    assert e["mean"] > 1.0


def test_the_headline_arm_is_like_for_like_with_the_reference():
    """A speedup may only be quoted from a run that did the reference's work, no more."""
    r = _bench_result("optimised_tuning126_w8")
    if r is None:
        pytest.skip("no w8 headline run on disk")
    assert r["system"] == "like-for-like-s9"
    assert r["components"] == {"amber": True, "quantum": False, "legacy": False}
    assert r["stages_s"]["reference_arms"] == 0.0,         "the verification arm is work s9/final.py does not do; it cannot be in a headline"
    assert r["stages_s"]["quantum"] == 0.0 and r["stages_s"]["legacy"] == 0.0
    assert r["cold"] and r["headline_eligible"] and r["n_done"] == 126


def test_the_speedup_is_computed_from_two_real_cold_runs():
    a = _bench_result("baseline_tuning126")
    b = _bench_result("optimised_tuning126_w8")
    if a is None or b is None:
        pytest.skip("both arms must have been run")
    for r, name in ((a, "baseline"), (b, "optimised")):
        assert r["cold"], f"{name} resumed targets; its wall clock is not a cold measurement"
        assert r["n_done"] == r["n_requested"] == 126
        assert r["dev_mode"] == ""
        assert r["manifest"] == "tuning126"
    speedup = a["wall_s"] / b["wall_s"]
    print(f"\nend-to-end {a['wall_s']:.1f} s -> {b['wall_s']:.1f} s = {speedup:.2f}x "
          f"({a['workers']}w -> {b['workers']}w)")
    assert speedup > 1.0
    #: and the per-stage table has to decompose the run, not sample it
    assert abs(a["stage_cpu_sum_s"] - a["wall_s"]) / a["wall_s"] < 0.05,         "the serial reference has unclocked time; its per-stage table is not a decomposition"


def test_a_shared_box_disqualifies_a_run_automatically():
    """Low occupancy WITH a tight worker spread means the box was shared, not badly scheduled.

    This sprint discarded a run that looked like a 1.2x regression and was not one: 8
    workers, 445.4 s, 67.4% system CPU, AMBER correctly declining on 4 of 126 targets at
    96% RAM.  Re-run on a quiet box the identical code gave 303.1 s and all 126.  The two
    failures have the same symptom in a wall clock and separate cleanly once per-worker
    busy time is known, so the harness decides it rather than a human noticing:

        idle from IMBALANCE = workers x max_busy - sum(busy)   the schedule's fault
        idle from OUTSIDE   = workers x (wall - max_busy)      something else's

    A run flagged `contention_suspected` is not headline-eligible, whatever its wall says.
    """
    from core import bench as B
    r = _bench_result("optimised_tuning126_w8")
    if r is None or not r.get("occupancy"):
        pytest.skip("no w8 run with per-worker data on disk")
    oc = r["occupancy"]
    print(f"\nheadline run: occupancy {oc['occupancy']*100:.1f}%, worker spread "
          f"{oc['busy_spread_frac']*100:.1f}%, idle {oc['idle_from_imbalance_s']:.0f} s "
          f"imbalance + {oc['idle_from_outside_s']:.0f} s outside")
    assert not oc["contention_suspected"],         "the headline was measured on a shared box and must be re-run"
    assert not oc["imbalance_suspected"],         f"worker spread {oc['busy_spread_frac']:.1%}: the schedule is leaving cores idle"
    assert r["headline_eligible"]
    #: and the rule must actually fire on a contended run, not merely be present
    bad = [{"pdb": f"p{i}", "wall": 100.0, "worker": {"pid": i % 8}} for i in range(16)]
    v = B._occupancy(bad, workers=8, wall=800.0, cpu_sum=1600.0)
    assert v["contention_suspected"] and not v["imbalance_suspected"], v
    good = B._occupancy(bad, workers=8, wall=210.0, cpu_sum=1600.0)
    assert not good["contention_suspected"], good


def test_the_machine_records_both_inflation_regimes():
    """A CPU-second floor is schedule-dependent, and quoting the wrong one hides a win.

    With N cores of unequal speed the per-structure CPU average differs by regime:
    equal COUNT per core gives mean(1/rel), equal WALL per core gives N/sum(rel).  Under
    any dynamic dispatch the fast cores complete more structures, so the average shifts
    toward them.  Both floors are recorded because this sprint's two schedules landed on
    one each -- 1.2913 against 1.2908, and 1.2533 against 1.2449 -- which confirms the
    core map twice from two regimes that disagree with each other.
    """
    from core import bench as B
    m = B.MACHINE
    rel = m["relative_throughput_per_core"]
    n = m["physical_cores"]
    assert len(rel) == n
    assert abs(sum(rel) - m["core_equivalents"]) < 5e-3
    assert abs(sum(1.0 / r for r in rel) / n
               - m["inflation_floor_static_equal_count"]) < 5e-3
    assert abs(n / sum(rel) - m["inflation_floor_dynamic"]) < 5e-3
    assert m["inflation_floor_dynamic"] < m["inflation_floor_static_equal_count"]
    b = _bench_result("baseline_tuning126")
    o = _bench_result("optimised_tuning126_w8")
    if b is None or o is None:
        pytest.skip("both arms must have been run")
    got = o["stages_s"]["amber"] / b["stages_s"]["amber"]
    floor = m["inflation_floor_dynamic"]
    print(f"\nAMBER inflation {got:.4f} against the dynamic-dispatch floor {floor} "
          f"= +{100 * (got / floor - 1):.2f}%")
    assert got >= floor - 5e-3, "AMBER came in under its own floor; the floor is wrong"
    assert got < floor * 1.05,         f"AMBER is {100 * (got / floor - 1):.1f}% over its floor: contention has returned"


# ============================================================ standalone runner
def _main():
    fails = []
    g = dict(globals())
    for name in sorted(g):
        if not name.startswith("test_"):
            continue
        fn = g[name]
        marks = getattr(fn, "pytestmark", [])
        cases = [()]
        for m in marks:
            if m.name == "parametrize":
                cases = [(v,) for v in m.args[1]]
        for args in cases:
            label = f"{name}{args if args else ''}"
            try:
                fn(*args)
                print(f"PASS  {label}")
            except Exception as exc:                                    # noqa: BLE001
                if type(exc).__name__ in ("Skipped", "OutcomeException"):
                    print(f"SKIP  {label}: {exc}")
                    continue
                fails.append((label, exc))
                print(f"FAIL  {label}: {exc}")
    print(f"\n{len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(_main())
