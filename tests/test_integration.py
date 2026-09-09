"""The four mandatory components are GENUINE, and the integrated path leaks nothing.

A scientific-integrity audit, not a formality. Each component is checked against an
INDEPENDENT reference wherever one can be built, because a component's own test agreeing
with the component proves only self-consistency:

  AMBER   ff14SB charge / sigma / epsilon and the GBn2 per-particle parameters compared
          against a System built here, from the same XML, on the same topology -- not
          against `core.amber`'s own numbers. Plus the pinned 1A13 interaction energy.
  Legacy  the 11-term decomposition compared column by column against the shipped
          `legacy_field`, with the weight vector's ORDER shown to be load-bearing and the
          term matrix shown to be non-degenerate.
  VQE     the ansatz's closed form compared against a `np.kron` statevector simulation --
          if the "circuit" is a classical shortcut, that is where it shows.
  CVaR    the tail, not the mean; and the corrected gradient, measured against
          finite differences computed from an independently rebuilt distribution.

Plus the standard leakage guard and the determinism/resumability checks.

The heavier checks (anything that builds an OpenMM System or runs a pipeline stage) are
gated on `VERIFY_SLOW=1`.
"""
from __future__ import annotations

import copy
import json
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "verify"))

import core                                                          # noqa: E402
from core import pipeline as P                                       # noqa: E402
from core import quantum as Q                                        # noqa: E402

slow = pytest.mark.skipif(os.environ.get("VERIFY_SLOW") != "1",
                          reason="set VERIFY_SLOW=1 for the OpenMM / pipeline checks")

#: 1A13's native interaction energy (nonbonded + solvation) at k=10, steps=0, tol=1.0.
REF_1A13_INTERACTION = -489.9138948277905
CANON_ALPHABET = "ARNDCQEGHILKMFPSTWYV"


# =========================================================== independent references
I2 = np.eye(2)


def _ry(t):
    c, s = np.cos(t / 2.0), np.sin(t / 2.0)
    return np.array([[c, -s], [s, c]])


def _op_on(n, q, g):
    m = np.array([[1.0]])
    for w in range(n):
        m = np.kron(m, g if w == q else I2)
    return m


def _cnot(n, ctrl, tgt):
    dim = 1 << n
    m = np.zeros((dim, dim))
    for j in range(dim):
        k = j ^ (1 << (n - 1 - tgt)) if (j >> (n - 1 - ctrl)) & 1 else j
        m[k, j] = 1.0
    return m


def _kron_probs(n, layers, theta, ring=True):
    """RY/CNOT circuit simulated with dense unitaries. Shares no code with `core`."""
    psi = np.zeros(1 << n)
    psi[0] = 1.0
    th = np.asarray(theta, float).reshape(layers, n)
    for L in range(layers):
        for q in range(n):
            psi = _op_on(n, q, _ry(th[L, q])) @ psi
        for q in range(n - 1):
            psi = _cnot(n, q, q + 1) @ psi
        if ring and n > 2:
            psi = _cnot(n, n - 1, 0) @ psi
    p = psi ** 2
    return p / p.sum()


def _loop_cvar(E, p, alpha):
    """CVaR of the lower alpha tail, accumulated by an explicit loop."""
    acc_m = acc_v = 0.0
    for i in sorted(range(len(E)), key=lambda i: (E[i], i)):
        take = min(p[i], alpha - acc_m)
        if take <= 0:
            break
        acc_v += take * E[i]
        acc_m += take
        if acc_m >= alpha:
            break
    return acc_v / alpha


def _cos(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d else float("nan")


# ================================================================= VQE genuineness
def test_vqe_ansatz_is_a_real_circuit_not_a_classical_shortcut():
    """`OneLayerAnsatz`'s closed form must BE the RY/CNOT circuit it claims to be.

    The ansatz replaces simulation with a Bernoulli-XOR product form. That is a legitimate
    optimisation only if it reproduces the circuit exactly; if it does not, the "quantum"
    distribution is a classical model wearing the name.
    """
    rng = np.random.default_rng(31337)
    for n in (4, 6, 8):
        for ring in (True, False):
            th = rng.uniform(0.2, 2.9, n)
            p_core = Q.OneLayerAnsatz(n, ring=ring).probs(th)
            p_ref = _kron_probs(n, 1, th, ring=ring and n > 2)
            assert np.max(np.abs(p_core - p_ref)) < 1e-12, (n, ring)


def test_statevector_circuit_matches_a_kron_simulation():
    """The batched statevector must equal a dense-unitary simulation of the same gates."""
    rng = np.random.default_rng(7)
    n, layers = 7, 3
    th = rng.normal(0, 0.8, n * layers)
    p_core = Q.StatevectorCircuit(n, layers, ring=True).probs(th)
    p_ref = _kron_probs(n, layers, th, ring=True)
    assert np.max(np.abs(p_core - p_ref)) < 1e-12


def test_the_state_is_parameterised_and_entangled():
    """A genuine variational state moves with theta and is not a product state."""
    rng = np.random.default_rng(5)
    circ = Q.StatevectorCircuit(7, 3)
    th = rng.normal(0, 0.7, circ.n_params())
    p0, p1 = circ.probs(th), circ.probs(th + 0.3)
    assert abs(p0.sum() - 1.0) < 1e-12
    assert np.max(np.abs(p1 - p0)) > 1e-3, "the distribution does not respond to theta"
    psi = np.sqrt(np.maximum(p0, 0)).reshape(1 << 3, 1 << 4)
    sv = np.linalg.svd(psi, compute_uv=False)
    assert sv[1] / sv[0] > 1e-6, "state factorises; the entangler is not doing anything"


# ================================================================ CVaR genuineness
def test_cvar_is_the_tail_not_the_mean():
    rng = np.random.default_rng(1)
    e = rng.standard_normal(4096)
    for a in (0.05, 0.15, 0.5):
        v = Q.cvar_from_samples(e, a)
        k = max(1, int(np.floor(a * e.size)))
        assert abs(v - float(np.sort(e)[:k].mean())) < 5e-3, a
        assert v < e.mean() - 0.5, (a, v, e.mean())


def test_parameter_shift_equals_independent_finite_differences():
    """The recorded signature of a correct gradient: cosine 1.000000.

    The finite differences here are taken on a CVaR computed by an explicit loop over a
    distribution rebuilt with `np.kron` -- no part of `core.quantum` is on the reference
    side, so this is a check and not a restatement.
    """
    rng = np.random.default_rng(20260904)
    n, layers, alpha = 6, 3, 0.15
    th = rng.normal(0, 0.8, n * layers)
    E = rng.normal(0, 1.0, 1 << n)
    circ = Q.StatevectorCircuit(n, layers)
    g_ps = Q.grad_cvar_paramshift(circ, th, E, alpha)

    h, g_ref = 1e-6, np.zeros(th.size)
    for k in range(th.size):
        tp, tm = th.copy(), th.copy()
        tp[k] += h
        tm[k] -= h
        g_ref[k] = (_loop_cvar(E, _kron_probs(n, layers, tp), alpha)
                    - _loop_cvar(E, _kron_probs(n, layers, tm), alpha)) / (2 * h)
    assert round(_cos(g_ps, g_ref), 6) == 1.000000, _cos(g_ps, g_ref)
    assert abs(np.linalg.norm(g_ps) / np.linalg.norm(g_ref) - 1.0) < 1e-6


def test_the_cvar_gradient_baseline_defect_has_not_returned():
    """THE regression test: `b(x) = mean * 1[x in tail]` is biased and must not be default.

    Measured with ZERO sampling noise -- the estimator is evaluated as an exact
    expectation over the enumerated register -- so any cosine below 1 is BIAS, with no
    variance left to blame. The constant baseline must land at 1.0; the tail-only form is
    computed alongside so the defect is priced rather than merely asserted absent.
    """
    rng = np.random.default_rng(20260904)
    n, layers, alpha = 7, 3, 0.15
    th = rng.normal(0, 0.8, n * layers)
    E = rng.normal(0, 1.0, 1 << n)
    circ = Q.StatevectorCircuit(n, layers)
    p = circ.probs(th)
    g_ps = Q.grad_cvar_paramshift(circ, th, E, alpha)
    q = Q.cvar_exact(E, p, alpha)[1]
    PR = circ.probs_batch(circ._shift_grid(th, np.pi / 2))
    G = ((PR[0::2] - PR[1::2]) / 2.0 / np.maximum(p, 1e-15)[None, :]).T

    w = np.where(E < q, (E - q) / alpha, 0.0)

    def expect(weights):
        return (p[:, None] * weights[:, None] * G).sum(0)

    g_const = expect(w - (p * w).sum())
    m = E < q
    wt = w.copy()
    wt[m] -= (p[m] * wt[m]).sum() / p[m].sum()
    wt[~m] = 0.0
    g_tail = expect(wt)

    assert round(_cos(g_const, g_ps), 6) == 1.000000, _cos(g_const, g_ps)
    assert abs(np.linalg.norm(g_const) / np.linalg.norm(g_ps) - 1.0) < 1e-6
    # the defect, measured: biased in direction AND magnitude
    assert _cos(g_tail, g_ps) < 0.80, _cos(g_tail, g_ps)
    assert np.linalg.norm(g_tail) / np.linalg.norm(g_ps) < 0.90


def test_the_shipped_default_baseline_is_the_constant_one():
    """A correct estimator that nothing calls is not a fix."""
    import inspect
    for name in ("cvar_gradient", "grad_cvar_score"):
        fn = getattr(Q, name)
        assert inspect.signature(fn).parameters["baseline"].default == "const", name
    rng = np.random.default_rng(2)
    an = Q.OneLayerAnsatz(12, ring=True)
    th = rng.uniform(0.3, 2.8, 12)
    bits = an.sample(th, 256, np.random.default_rng(9))
    E = rng.standard_normal(256)
    assert np.array_equal(Q.cvar_gradient(an, th, bits, E, 0.25)[0],
                          Q.cvar_gradient(an, th, bits, E, 0.25, baseline="const")[0])


def test_the_deployed_vqe_path_uses_the_exact_gradient():
    """`run_cvar_vqe` optimises the parameter-shift gradient, which HAS no baseline.

    That is what makes the historical defect structurally unreachable on the deployed
    path, so it is asserted rather than assumed.
    """
    import inspect
    src = inspect.getsource(Q.free_energy)
    assert "_shift_grid" in src and "baseline" not in src
    assert "free_energy" in inspect.getsource(Q.run_cvar_vqe)
    E = np.random.default_rng(0).standard_normal(128)
    p, cv, H, circ = Q.run_cvar_vqe(E, alpha=0.25, T=0.3, n=7, layers=3, iters=5, seed=0)
    assert abs(float(np.sum(p)) - 1.0) < 1e-9
    assert np.isfinite(cv) and np.isfinite(H) and H > 0.0


# ============================================================== Legacy genuineness
def test_legacy_is_the_real_eleven_term_field_with_its_fitted_weights():
    import energy_terms as et
    lf = core.backend("legacy")
    assert list(lf.TERMS) == list(et.TERM_NAMES)
    assert len(lf.TERMS) == 11
    import legacy_field as ref
    assert ({k: float(v) for k, v in lf.FITTED_WEIGHTS.items()}
            == {k: float(v) for k, v in ref.FITTED_WEIGHTS.items()})
    assert len(set(lf.FITTED_WEIGHTS.values())) > 1, "weights are a constant, not fitted"
    assert set(lf.FITTED_WEIGHTS) == set(lf.TERMS)


@slow
def test_legacy_terms_match_the_shipped_module_on_a_real_pool():
    """Column by column on real retrieved geometry, and non-degenerate."""
    import legacy_field as ref
    import protein_geometry as geo
    from verify.legacy_audit import RepShim, real_pool          # noqa: F401
    lf = core.backend("legacy")
    p, PHI, PSI = real_pool()
    BB = geo.build_backbone_batch(PHI, PSI)
    T_new = np.asarray(lf.BatchLegacy(p.seq, RepShim(PHI, PSI))
                       .terms_from_coords(BB, phi=PHI, psi=PSI), float)
    T_ref = np.asarray(ref.BatchLegacy(p.seq, RepShim(PHI, PSI))
                       .terms_from_coords(BB, phi=PHI, psi=PSI), float)
    assert np.array_equal(T_new, T_ref)
    live = T_new.std(0) > 1e-12
    Z = (T_new[:, live] - T_new[:, live].mean(0)) / T_new[:, live].std(0)
    assert np.linalg.matrix_rank(Z, tol=1e-8) >= int(live.sum()) - 1


# =============================================================== AMBER genuineness
@slow
def test_amber_is_real_ff14sb_gbn2_with_the_shipped_parameters():
    """Not just "a CustomGBForce is present" -- the PARAMETERS must be ff14SB's.

    A System can carry a real `CustomGBForce` whose per-particle parameters have been
    replaced. So charge, sigma, epsilon and the GBn2 parameters are compared against a
    System built here from the same XML on the same topology.
    """
    from openmm import app, unit
    import core.amber as A
    from verify.amber_audit import build_native

    p, rep, _ = build_native()
    H = A.builder_for(p.seq, rep)
    names = sorted(f.__class__.__name__ for f in H.system.getForces())
    assert "CustomGBForce" in names and "NonbondedForce" in names
    assert "amber14/protein.ff14SB.xml" in H.describe()["forcefield"]
    assert "implicit/gbn2.xml" in H.describe()["forcefield"]

    ref_sys = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml") \
        .createSystem(H.topology, nonbondedMethod=app.NoCutoff, constraints=None,
                      rigidWater=False, implicitSolventKappa=0.0 / unit.nanometer)
    nb = next(f for f in H.system.getForces()
              if f.__class__.__name__ == "NonbondedForce")
    nbr = next(f for f in ref_sys.getForces()
               if f.__class__.__name__ == "NonbondedForce")
    for i in range(nb.getNumParticles()):
        a, b = nb.getParticleParameters(i), nbr.getParticleParameters(i)
        for x, y in zip(a, b):
            assert x == y, (i, x, y)

    gb = next(f for f in H.system.getForces()
              if f.__class__.__name__ == "CustomGBForce")
    gbr = next(f for f in ref_sys.getForces()
               if f.__class__.__name__ == "CustomGBForce")
    assert gb.getNumComputedValues() == 2 and gb.getNumPerParticleParameters() == 7
    for i in range(gb.getNumParticles()):
        assert list(gb.getParticleParameters(i)) == list(gbr.getParticleParameters(i)), i


@slow
def test_the_1a13_native_interaction_energy_is_bit_exact():
    import core.amber as A
    from verify.amber_audit import build_native
    p, rep, rb = build_native()
    cm = A.refine_coords(p.seq, rep, rb, k_restraint=A.K_MODERATE, steps=0,
                         tolerance=1.0, components=True)["components"]
    assert cm["nonbonded"] + cm["solvation"] == REF_1A13_INTERACTION


@slow
def test_single_point_energy_is_invariant_under_rigid_translation():
    """Real physics has no absolute position; a hash of the coordinates would.

    The bar is RELATIVE, and deliberately so. The unminimised native carries ~2e8
    kcal/mol of builder strain, and differencing two such numbers loses ~8 digits to
    cancellation before any physics is involved; an absolute tolerance here would be
    measuring float64's mantissa, not translation invariance. 1e-5 relative is ~3 orders
    tighter than any real position dependence would be and still above the noise.
    """
    import core.amber as A
    from verify.amber_audit import build_native
    p, rep, rb = build_native()
    kw = dict(k_restraint=0.0, steps=-1, tolerance=1e9, components=True)
    a = A.refine_coords(p.seq, rep, rb, **kw)
    b = A.refine_coords(p.seq, rep,
                        {k: v + np.array([1.75, -0.5, 0.25]) for k, v in rb.items()}, **kw)
    ea, eb = float(a["energy"]), float(b["energy"])
    rel = abs(ea - eb) / max(abs(ea), 1.0)
    print(f"\nsingle point: {ea!r} vs {eb!r}  relative {rel:.3e}")
    assert rel < 1e-5, (ea, eb, rel)
    # the interaction terms are the ones the science reads, and they are small enough
    # that an absolute bar is meaningful on them
    for term in ("nonbonded", "solvation"):
        x, y = a["components"][term], b["components"][term]
        assert abs(x - y) / max(abs(x), 1.0) < 1e-4, (term, x, y)


# ==================================================================== hazards
def test_only_the_canonical_alphabet_is_used_to_encode():
    """`ARNDCQEGHILKMFPSTWYV` is canonical; the other order must never encode anything."""
    d = core.backend("data")
    assert d.ALPHABET == CANON_ALPHABET
    assert sorted(d.ALPHABET_ALT) == sorted(CANON_ALPHABET)
    assert d.ALPHABET_ALT != d.ALPHABET
    # the encoder must agree with the canonical order, position for position
    for i, a in enumerate(CANON_ALPHABET):
        assert int(d.encode(a)[0]) == i, (a, i)


def test_identity_still_normalises_by_the_longer_sequence():
    """Leaky BY DESIGN. "Fixing" it would silently move the benchmark's target set."""
    d = core.backend("data")
    assert abs(d.identity("ACDEFG", "ACDEFGHIKL") - 0.6) < 1e-12
    assert d.identity("ACDEFG", "ACDEFGHIKL") == d.identity("ACDEFGHIKL", "ACDEFG")
    assert d.IDENTITY_THRESHOLD == 0.6


def test_the_ranking_paths_use_a_stable_argsort():
    """Unstable sorting moves pool membership by up to 47 of 500 on tied similarities."""
    assert P.PROD.tie_break == "stable"
    import inspect
    for fn in (P.retrieve, P.filter_pool):
        src = inspect.getsource(fn)
        assert "argsort" not in src or "cfg.tie_break" in src, fn.__name__
    sc = np.round(np.random.default_rng(3).standard_normal(5000), 2)
    a = np.argsort(sc, kind="stable")[:500]
    b = np.argsort(sc, kind="stable")[:500]
    assert np.array_equal(a, b)


def test_batched_contraction_is_alignment_invariant():
    """BLAS picks its SIMD path on buffer alignment; the result must not depend on it."""
    rng = np.random.default_rng(11)
    A0 = rng.standard_normal((400, 96))
    v = rng.standard_normal(96)
    ref = A0 @ v
    for off in range(1, 9):
        buf = np.empty((400, 96 + off))
        buf[:, off:] = A0
        assert np.array_equal(buf[:, off:] @ v, ref), off


# ==================================================================== leakage
@slow
@pytest.mark.parametrize("pdb", ["1CS9", "1CB3"])
def test_nan_poisoning_the_native_changes_no_deployable_quantity(pdb):
    """THE leakage guard. If a native reaches a deployable stage, the NaN propagates."""
    db = core.backend("data")
    folds = db.folds(P.PROD.n_folds)
    t = [x for x in P.manifest("smoke8") if x.pdb == pdb][0]
    fold = int(folds[t.seq])

    clean, _, _ = P.run_target(t, fold, P.PROD)
    dirty_t = copy.deepcopy(t)
    for f in ("ca", "phi", "psi"):
        if hasattr(dirty_t, f):
            object.__setattr__(dirty_t, f,
                               np.full(np.asarray(getattr(t, f), float).shape, np.nan))
    dirty, _, _ = P.run_target(dirty_t, fold, P.PROD)

    def flat(rec):
        out = {}
        for k, v in rec.items():
            if isinstance(v, np.ndarray):
                out[k] = np.asarray(v, float)
            elif isinstance(v, dict):
                for k2, v2 in v.items():
                    if isinstance(v2, np.ndarray):
                        out[f"{k}.{k2}"] = np.asarray(v2, float)
                    elif isinstance(v2, (int, float)):
                        out[f"{k}.{k2}"] = np.asarray([v2], float)
            elif isinstance(v, (int, float)) and not isinstance(v, bool):
                out[k] = np.asarray([v], float)
        return out

    a, b = flat(clean), flat(dirty)
    keys = sorted(set(a) & set(b))
    assert len(keys) >= 12, keys
    for k in keys:
        assert not np.isnan(b[k]).any(), f"NaN leaked into {k}"
        assert np.array_equal(a[k], b[k]), f"{k} changed when the native was poisoned"


# ============================================================ determinism / cache
@slow
def test_the_same_target_is_bit_identical_across_processes():
    from verify.determinism_audit import run_child
    a, b = run_child("1CS9"), run_child("1CS9")
    keys = sorted(set(a) & set(b))
    assert len(keys) >= 12
    assert all(a[k] == b[k] for k in keys), [k for k in keys if a[k] != b[k]]


@slow
def test_the_result_does_not_depend_on_the_ambient_thread_count():
    """The hazard that moved single-point energies by 144 kcal/mol."""
    from verify.determinism_audit import run_child
    a = run_child("1CS9")
    c = run_child("1CS9", {"OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "4",
                           "MKL_NUM_THREADS": "4", "OPENMM_CPU_THREADS": "4"})
    keys = sorted(set(a) & set(c))
    assert all(a[k] == c[k] for k in keys), [k for k in keys if a[k] != c[k]]


def test_cache_keys_discriminate_every_type_json_would_conflate():
    from core import cache as C
    ns = "verify_types"
    ks = [C.key(ns, 1, x=1), C.key(ns, 1, x=1.0), C.key(ns, 1, x=True),
          C.key(ns, 1, x="1"), C.key(ns, 1, x=None), C.key(ns, 1, x=[1, 2]),
          C.key(ns, 1, x=(1, 2)), C.key(ns, 1, x=np.nextafter(1.0, 2.0)),
          C.key(ns, 2, x=1), C.key(ns, 1, x=1, y=0)]
    assert len(set(ks)) == len(ks), "two distinct parameter sets share a cache key"
    with pytest.raises(TypeError):
        C.key(ns, 1, x=object())


def test_a_key_missing_a_parameter_raises_instead_of_serving_a_stale_array():
    """Forge the sidecar -- the signature of an incomplete key -- and demand a raise."""
    from core import cache as C
    ns = "verify_collide_test"
    k = C.store(ns, 1, {"a": np.arange(4.0)}, alpha=1, beta="x")
    _, side = C._paths(ns, k)
    C._MEM.clear()
    meta = json.load(open(side))
    meta["params"] = json.dumps(C._canon({"alpha": 2, "beta": "x"}),
                                sort_keys=True, separators=(",", ":"))
    json.dump(meta, open(side, "w"))
    C._MEM.clear()
    try:
        with pytest.raises(C.CacheCollision):
            C.load(ns, 1, alpha=1, beta="x")
    finally:
        C.clear(ns)


def test_the_projection_gradient_is_in_the_cache_key():
    """PROJECT_GRAD changes the emitted structure, so it must change the cache key.

    `core/project.py` reads its gradient from an environment variable at import time:

        GRAD = os.environ.get("PROJECT_GRAD", "analytic")

    and its own docstring says the two settings "land on a different point of the same
    basin". Measured directly on smoke8 targets with the pipeline cache bypassed, `fd` and
    `analytic` produce CA traces up to 0.45 A apart and post-AMBER traces up to 0.45 A
    apart. But `Config` has no field for it, so `Config.key()` cannot see it and both
    settings address the same per-target cache directory.

    Observed consequence: an `fd` run of smoke8 computed 8 targets fresh; the `analytic`
    run that followed reported the SAME cfg_key and finished all 8 in 0.215 s, serving the
    `fd` arm's arrays under the analytic arm's name. That is a stale-result collision of
    exactly the kind `core/cache.py` raises on -- but the per-target pipeline cache is a
    separate mechanism keyed only on `Config.key()`, with no parameter sidecar to check.

    The fix is to put the gradient selection in `Config` (a field, hashed like every other
    science parameter) rather than in a module-level environment global.
    """
    import core.project as cp
    keyed = any("grad" in f.lower() for f in P.PROD.science())
    assert keyed, (
        "PROJECT_GRAD selects between two projections that emit different structures "
        f"(module default {cp.GRAD!r}) but no Config field carries it, so two runs whose "
        "results differ share a cache key and the second serves the first's arrays. "
        f"Config science fields: {sorted(P.PROD.science())}")


def test_every_science_parameter_moves_the_config_key():
    """A parameter that changes a number but not the key would serve a stale result."""
    from dataclasses import replace, fields
    base = P.PROD.key({"x": 1})
    science = set(P.PROD.science())
    for f in fields(P.PROD):
        if f.name not in science:
            continue
        cur = getattr(P.PROD, f.name)
        if isinstance(cur, bool):
            new = not cur
        elif isinstance(cur, int):
            new = cur + 1
        elif isinstance(cur, float):
            new = cur + 1.0
        elif isinstance(cur, str):
            new = cur + "_x"
        else:
            continue
        assert replace(P.PROD, **{f.name: new}).key({"x": 1}) != base, f.name
    # and the backend set is part of the key, so an optimised result can never be
    # served out of a baseline cache
    assert P.PROD.key({"x": 1}) != P.PROD.key({"x": 2})
