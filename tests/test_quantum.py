"""Equivalence and genuineness tests for `core.quantum`.

What this suite has to support
Two claims, and neither of them is "the consolidation looks right".

**1. It is still a real VQE with a real CVaR.**  A consolidation of a quantum component
has one catastrophic failure mode -- quietly replacing the circuit with classical
enumeration or the CVaR with a mean -- and it would not show up in a speed number or in a
diff review.  So the ansatz families are checked against an INDEPENDENT dense statevector
built here from Kronecker products (`_dense_state`), which shares no code with
`core.quantum`; the CVaR is checked to be the tail and not the mean; and the search path
is checked to touch a number of bitstrings that is a vanishing fraction of ``2**n``.

**2. The floats did not move.**  Every routine that was rewritten is compared against the
module it was lifted from -- `qansatz`, `foldvqe`, `s8/integrate` -- on the same inputs
with the same seeds.  Where the answer is bit-identical the assertion is ``==``, not
``approx``, because that is a stronger and more honest test than a tolerance.  Where the
reduction ORDER genuinely changed (the MPS contraction is now log-depth) the difference is
MEASURED and pinned with a stated tolerance, and the tolerance is justified in the test.

The one test that matters most
`test_cvar_gradient_baseline_is_constant_not_tail_only` pins a defect that was already
found and fixed once: `qansatz.cvar_gradient` subtracts the tail mean from the tail entries
only, i.e. it uses ``b(x) = m * 1[x in tail]``, which is a FUNCTION of x and therefore not
an unbiased control variate.  That test does not read the source to check which baseline is
in use -- it MEASURES both against an exact classical reference on an enumerable register
and asserts the shipped default lands on +1.000000 rather than +0.655634.  Reintroducing
the defect makes it fail with the number that names the bug.

Run:  python -m pytest tests/test_quantum.py -q
      python tests/test_quantum.py          (no-pytest fallback runner)
"""
import importlib.util
import math
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import numpy as np                                                    # noqa: E402

import core.quantum as Q                                              # noqa: E402
import foldvqe as FV                                                  # noqa: E402
import qansatz                                                        # noqa: E402
import vqe as VQE                                                     # noqa: E402

Q.limit_threads(2)

#: The s9 `stage_grad` configuration, reproduced verbatim so the recorded numbers are
#: reproduced rather than approximated. SEED and the draw order come from `s9/refine.py`.
S9_SEED = 20260914
S9_NQ = 10
S9_TRIALS = 12
S9_SHOTS = 4096
S9_ALPHAS = (0.10, 0.25, 0.50)

#: The recorded audit numbers this file exists to keep reproducing.
REC_CONST_COS = 1.000000
REC_TAIL_COS = 0.655634
REC_TAIL_NORM = 0.758


def _s8():
    """`s8/integrate.py` loaded by path -- it is not an importable package member."""
    spec = importlib.util.spec_from_file_location(
        "_s8_integrate", os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "s8", "integrate.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _cos(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-300))


# ============================================================ independent reference
def _dense_state(theta, n, nblocks, layers, ring, entangler="cnot"):
    """A dense ``2**n`` statevector built from Kronecker products.

    Deliberately shares NO code with `core.quantum`: gates are assembled as full
    ``2**n x 2**n`` matrices and multiplied, which is the least clever and least
    coincidence-prone way to say what the circuit means.  It is the ground truth that both
    the MPS and the batched statevector are checked against, so "the ansatz is still the
    circuit it claims to be" is a measured statement here rather than an assumption.
    """
    th = np.asarray(theta, float).reshape(nblocks, n)
    dim = 1 << n
    psi = np.zeros(dim)
    psi[0] = 1.0

    def kron_ry(q, t):
        c, s = math.cos(t / 2.0), math.sin(t / 2.0)
        g = np.array([[c, -s], [s, c]])
        M = np.eye(1)
        for w in range(n):
            M = np.kron(M, g if w == q else np.eye(2))
        return M

    def cnot(ctrl, tgt):
        M = np.zeros((dim, dim))
        for i in range(dim):
            cb = (i >> (n - 1 - ctrl)) & 1
            j = i ^ (1 << (n - 1 - tgt)) if cb else i
            M[j, i] = 1.0
        return M

    for l in range(nblocks):
        for q in range(n):
            psi = kron_ry(q, th[l, q]) @ psi
        if l >= layers or entangler == "none":
            continue
        for q in range(n - 1):
            psi = cnot(q, q + 1) @ psi
        if ring and n > 2:
            psi = cnot(n - 1, 0) @ psi
    return psi


def _brute_cvar(energies, probs, alpha):
    """CVaR by walking the sorted mass ladder in Python. The definition, transcribed."""
    e = np.asarray(energies, float)
    p = np.asarray(probs, float) / np.sum(probs)
    order = np.argsort(e, kind="stable")
    acc, tot = 0.0, 0.0
    for k in order:
        take = min(p[k], alpha - acc)
        if take <= 0:
            break
        tot += take * e[k]
        acc += take
        if acc >= alpha:
            break
    return tot / acc


# =============================================== 1. exact-enumeration agreement
def test_cvar_objective_matches_exact_enumeration():
    """The sampled CVaR objective converges to the exact CVaR of the circuit's own
    distribution, and that distribution is the circuit's.

    Three things are checked on ONE small register, because they are three links of one
    chain and checking them separately would not prove the chain:

      1. the ansatz's ``probs`` equals ``|psi|^2`` of the independent dense simulation
         (the state is the circuit's);
      2. ``cvar_from_probs`` on that distribution equals a brute-force walk of the sorted
         mass ladder (the functional is CVaR);
      3. the sampled ``cvar_from_samples`` over shots drawn from the ansatz converges to
         it (what the optimiser actually evaluates is an estimate of that quantity).
    """
    n = 8
    rng = np.random.default_rng(11)
    E = rng.standard_normal(1 << n)
    for kind, ans, nb, lay, ring in (
            ("chain", Q.OneLayerAnsatz(n, ring=True), 1, 1, True),
            ("chain_ry", Q.MPSAnsatz(n, 1, final_ry=True), 2, 1, False),
            ("deep", Q.MPSAnsatz(n, 2, final_ry=True), 3, 2, False)):
        theta = rng.uniform(0.3, math.pi - 0.3, ans.n_params())
        p_ans = ans.probs(theta)
        p_ref = _dense_state(theta, n, nb, lay, ring) ** 2
        assert np.abs(p_ans - p_ref).max() < 1e-12, kind
        assert abs(p_ans.sum() - 1.0) < 1e-12, kind

        for alpha in (0.05, 0.25, 1.0):
            exact = Q.cvar_from_probs(E, p_ref, alpha)[0]
            assert abs(exact - _brute_cvar(E, p_ref, alpha)) < 1e-12, (kind, alpha)
            # and the mean is a DIFFERENT number, so "CVaR" is not a misnamed mean
            if alpha < 1.0:
                assert exact < float(p_ref @ E) - 1e-6, (kind, alpha)

            bits = ans.sample(theta, 400_000, np.random.default_rng(5))
            idx = bits.astype(np.int64) @ (1 << np.arange(n - 1, -1, -1))
            samp = Q.cvar_from_samples(E[idx], alpha)
            assert abs(samp - exact) < 0.05, (kind, alpha, samp, exact)


def test_cvar_is_the_tail_and_the_boundary_cases_are_right():
    """alpha=1 is the mean; alpha -> 0 is the minimum; in between it is strictly ordered."""
    rng = np.random.default_rng(3)
    e = rng.standard_normal(1000)
    assert abs(Q.cvar_from_samples(e, 1.0) - e.mean()) < 1e-12
    assert Q.cvar_from_samples(e, 1e-6) == e.min()
    vals = [Q.cvar_from_samples(e, a) for a in (0.01, 0.05, 0.2, 0.5, 1.0)]
    assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1)), vals


# ======================================== 2. parameter shift vs finite differences
def test_parameter_shift_matches_exact_finite_differences():
    """Reproduce the recorded ``cos 1.000000``.

    ``p(x) = <psi|Pi_x|psi>`` is a projector expectation, so each basis probability obeys
    the exact two-term shift rule for an RY generator; chained with ``dCVaR/dp`` that gives
    the analytic gradient.  Central finite differences on the exact CVaR value use none of
    that machinery -- no projector algebra, no envelope theorem -- so agreement between
    them is a real check rather than a restatement.
    """
    rng = np.random.default_rng(0)
    for n, layers in ((7, 3), (9, 2)):
        circ = Q.StatevectorCircuit(n, layers)
        E = rng.standard_normal(circ.dim)
        for alpha in (0.10, 0.25, 0.50):
            th = rng.normal(0, 1.0, circ.n_params())
            gp = Q.grad_cvar_paramshift(circ, th, E, alpha)
            gf = Q.grad_cvar_fd(circ, th, E, alpha)
            assert round(_cos(gp, gf), 6) == 1.000000, (n, alpha, _cos(gp, gf))
            rel = np.linalg.norm(gp - gf) / max(np.linalg.norm(gp), 1e-300)
            assert rel < 1e-7, (n, alpha, rel)


def test_finite_difference_reference_is_independent_machinery():
    """The FD reference must not be secretly calling the shift rule.

    A reference that reuses the routine it audits proves nothing.  This asserts the
    numerical fingerprint of that independence: FD carries a truncation error that scales
    with the step, and the shift rule does not.  If `grad_cvar_fd` were delegating, halving
    ``h`` would change nothing.
    """
    rng = np.random.default_rng(4)
    circ = Q.StatevectorCircuit(7, 3)
    E = rng.standard_normal(circ.dim)
    th = rng.normal(0, 1.0, circ.n_params())
    gp = Q.grad_cvar_paramshift(circ, th, E, 0.25)
    e_coarse = np.linalg.norm(Q.grad_cvar_fd(circ, th, E, 0.25, h=1e-3) - gp)
    e_fine = np.linalg.norm(Q.grad_cvar_fd(circ, th, E, 0.25, h=1e-5) - gp)
    assert e_coarse > 10 * e_fine, (e_coarse, e_fine)


# ============================== 3. The baseline regression test (most important)
def _stage_grad(nq=S9_NQ, trials=S9_TRIALS, shots=S9_SHOTS, alphas=S9_ALPHAS):
    """The s9 `stage_grad` audit, reproduced against an EXACT classical reference.

    At ``nq`` qubits all ``2**nq`` amplitudes of ``p_theta`` are closed form, so the CVaR of
    the exact distribution and its central-finite-difference gradient are exact.  Four
    estimators of the same quantity are compared to that reference:

      * EXACT-expectation score function with a CONSTANT baseline  (the fix),
      * EXACT-expectation score function with the shipped TAIL-ONLY baseline (the defect),
      * SAMPLED score function with each of the two.

    The exact-expectation rows are the point: they carry NO sampling noise, so a
    disagreement there is bias and nothing else.
    """
    rng = np.random.default_rng(S9_SEED + 11)
    an = Q.OneLayerAnsatz(nq, ring=True)
    allbits = Q.all_bitstrings(nq)

    def exact_p(theta):
        return np.exp(an.logp(theta, allbits))

    def wcvar(E, p, alpha):
        o = np.argsort(E, kind="stable")
        cw = np.cumsum(p[o])
        k = min(int(np.searchsorted(cw, alpha)) + 1, len(o))
        tail = np.zeros(len(E), bool)
        tail[o[:k]] = True
        return float(E[o[:k]].mean()), float(E[o[k - 1]]), tail

    def exact_cvar(theta, E, alpha):
        p = exact_p(theta)
        o = np.argsort(E, kind="stable")
        cw = np.cumsum(p[o])
        k = min(int(np.searchsorted(cw, alpha)) + 1, len(o))
        w = p[o[:k]].copy()
        w[-1] -= (cw[k - 1] - alpha)
        w = np.maximum(w, 0.0)
        return float((w * E[o[:k]]).sum() / max(w.sum(), 1e-300))

    out = {k: [] for k in ("exact_const", "exact_tail", "samp_const", "samp_tail")}
    norms = {k: [] for k in out}
    for _t in range(trials):
        for alpha in alphas:
            theta = rng.uniform(0.3, math.pi - 0.3, nq)
            E = rng.standard_normal(1 << nq)
            h = 1e-5
            gref = np.zeros(nq)
            for q in range(nq):
                tp, tm = theta.copy(), theta.copy()
                tp[q] += h
                tm[q] -= h
                gref[q] = (exact_cvar(tp, E, alpha) - exact_cvar(tm, E, alpha)) / (2 * h)
            p = exact_p(theta)
            G = an.grad_logp(theta, allbits)
            idx = rng.choice(1 << nq, size=shots, p=p / p.sum())
            sb, se = allbits[idx], E[idx]

            def sf(f):
                return (p[:, None] * f[:, None] * G).sum(0)

            _v, q0, tail = wcvar(E, p, alpha)
            f = np.zeros(len(E))
            f[tail] = -(q0 - E[tail]) / alpha
            g = {"exact_const": sf(f - (p * f).sum())}
            f2 = f.copy()
            f2[tail] -= (p[tail] * f[tail]).sum() / max(p[tail].sum(), 1e-300)
            f2[~tail] = 0.0
            g["exact_tail"] = sf(f2)
            g["samp_const"] = Q.cvar_gradient(an, theta, sb, se, alpha, "const")[0]
            g["samp_tail"] = Q.cvar_gradient(an, theta, sb, se, alpha, "tail")[0]
            for k, v in g.items():
                out[k].append(_cos(v, gref))
                norms[k].append(np.linalg.norm(v) / max(np.linalg.norm(gref), 1e-300))
    return ({k: float(np.mean(v)) for k, v in out.items()},
            {k: float(np.mean(v)) for k, v in norms.items()})


def test_cvar_gradient_baseline_is_constant_not_tail_only():
    """The regression test for the recorded CVaR-gradient baseline defect.

    A baseline ``b`` leaves a score-function estimator unbiased iff it is CONSTANT in x,
    because the correction term is ``b * E_p[grad log p] = 0``.  The shipped
    `qansatz.cvar_gradient` subtracts the tail mean from the tail entries and leaves the
    rest at zero -- ``b(x) = m * 1[x in tail]``, a function of x -- and the identity does
    not hold on a data-dependent subset.

    This test does not inspect the source to find out which baseline is in use.  It
    MEASURES both against an exact reference on an enumerable register and asserts that the
    shipped default lands on the fix's number.  If the defect is ever reintroduced, the
    assertion fails carrying +0.655634 -- the number that names the bug.
    """
    cos, nr = _stage_grad()
    # the fix: exact to round-off
    assert round(cos["exact_const"], 6) == REC_CONST_COS, cos
    assert abs(nr["exact_const"] - 1.0) < 1e-4, nr
    # the defect, reproduced and still priced where the record says it is
    assert abs(cos["exact_tail"] - REC_TAIL_COS) < 2e-3, cos
    assert abs(nr["exact_tail"] - REC_TAIL_NORM) < 2e-3, nr
    # and it is BIAS, not variance: the two rows above used the exact expectation, so
    # they contain no sampling noise at all. The sampled rows only confirm the ordering.
    assert cos["samp_const"] > 0.97, cos
    assert cos["samp_tail"] < 0.80, cos
    assert cos["samp_const"] - cos["samp_tail"] > 0.25, cos


def test_shipped_default_baseline_is_the_constant_one():
    """`cvar_gradient` called the way the search calls it must BE the constant form.

    The previous test proves the constant baseline is correct; this one proves it is what
    you get by default.  Both are needed: a correct estimator that nothing calls is not a
    fix.
    """
    rng = np.random.default_rng(2)
    for ans, npar in ((Q.OneLayerAnsatz(14, ring=True), 14),
                      (Q.MPSAnsatz(12, 1, final_ry=True), 24)):
        theta = rng.uniform(0.3, 2.8, npar)
        bits = ans.sample(theta, 256, np.random.default_rng(9))
        E = rng.standard_normal(256)
        g_default = Q.cvar_gradient(ans, theta, bits, E, 0.25)[0]
        g_const = Q.cvar_gradient(ans, theta, bits, E, 0.25, baseline="const")[0]
        g_tail = Q.cvar_gradient(ans, theta, bits, E, 0.25, baseline="tail")[0]
        assert np.array_equal(g_default, g_const)
        assert not np.allclose(g_default, g_tail, rtol=1e-3, atol=1e-9)


def test_tail_baseline_reproduces_the_shipped_routine_exactly():
    """``baseline="tail"`` must reproduce `qansatz.cvar_gradient` to round-off.

    This is what makes the defect MEASURABLE rather than merely described: the audit prices
    the code that actually shipped, not a paraphrase of it.  The recorded standard is
    5.6e-17 absolute; the assertion below is the relative form of it.
    """
    rng = np.random.default_rng(1)
    worst = 0.0
    for nq, layers in ((16, 1), (14, 2), (10, 1)):
        for cls in ("one", "mps"):
            if cls == "one":
                a_old = qansatz.OneLayerAnsatz(nq, ring=True)
                a_new = Q.OneLayerAnsatz(nq, ring=True)
            else:
                a_old = qansatz.MPSAnsatz(nq, layers, final_ry=True)
                a_new = Q.MPSAnsatz(nq, layers, final_ry=True)
            theta = rng.normal(0.0, 0.8, a_old.n_params()) + math.pi / 2
            bits = a_old.sample(theta, 256, np.random.default_rng(4))
            E = rng.standard_normal(256)
            for alpha in (0.10, 0.25, 0.50):
                g0, v0 = qansatz.cvar_gradient(a_old, theta, bits, E, alpha)
                g1, v1 = Q.cvar_gradient(a_new, theta, bits, E, alpha, baseline="tail")
                assert v0 == v1, (cls, nq, alpha)
                rel = np.abs(g0 - g1).max() / max(np.abs(g0).max(), 1e-300)
                worst = max(worst, rel)
    assert worst < 1e-14, worst


def test_constant_baseline_does_not_change_the_expectation():
    """A control variate must be exactly that: variance down, expectation unmoved.

    ``baseline="none"`` and ``baseline="const"`` differ by ``b * mean_i grad log p(x_i)``,
    which converges to zero.  Over many independent batches their means must agree, while
    the constant-baseline estimator's spread must be no worse.  The tail form fails the
    first half of that, which is the whole point.
    """
    rng = np.random.default_rng(7)
    nq = 9
    ans = Q.OneLayerAnsatz(nq, ring=True)
    theta = rng.uniform(0.4, 2.7, nq)
    E = rng.standard_normal(1 << nq)
    allbits = Q.all_bitstrings(nq)
    g_exact = Q.cvar_gradient_exact(ans, theta, allbits, E, 0.2)[0]
    draw = np.random.default_rng(23)
    acc = {"none": [], "const": [], "tail": []}
    for _ in range(300):
        bits = ans.sample(theta, 256, draw)
        e = E[bits.astype(np.int64) @ (1 << np.arange(nq - 1, -1, -1))]
        for b in acc:
            acc[b].append(Q.cvar_gradient(ans, theta, bits, e, 0.2, baseline=b)[0])
    m = {b: np.mean(v, axis=0) for b, v in acc.items()}
    assert _cos(m["none"], g_exact) > 0.95, _cos(m["none"], g_exact)
    assert _cos(m["const"], g_exact) > 0.95, _cos(m["const"], g_exact)
    # the tail form's MEAN over 300 independent batches is still far off: bias, not noise
    assert _cos(m["tail"], g_exact) < _cos(m["const"], g_exact) - 0.1, \
        (_cos(m["tail"], g_exact), _cos(m["const"], g_exact))


# ======================================== 5. sort vs partition tail equivalence
def test_sort_and_partition_select_the_same_tail():
    """``np.argpartition`` selection == stable-``argsort`` prefix, exactly.

    CVaR needs a SELECTION, not an ordering.  A stable sort places equal elements in index
    order, so ``argsort(e)[:k]`` is ``{e < q}`` plus the lowest INDICES among the ties at
    ``q`` -- which is what `tail_indices` constructs in O(N).  Asserted on the mask, the
    quantile and the value, with ``==`` rather than a tolerance.
    """
    rng = np.random.default_rng(5)
    # sizes deliberately straddle CVAR_SORT_CUTOFF, because both branches ship
    assert Q.CVAR_SORT_CUTOFF > 2
    for N in (1, 2, 7, Q.CVAR_SORT_CUTOFF - 1, Q.CVAR_SORT_CUTOFF,
              Q.CVAR_SORT_CUTOFF + 1, 4096):
        for alpha in (1.0, 0.5, 0.25, 0.1, 0.01, 1e-9):
            for trial in range(20):
                e = rng.standard_normal(N)
                old = qansatz.cvar(e, alpha)
                new = Q.cvar(e, alpha)
                assert old[0] == new[0], (N, alpha, trial)
                assert old[1] == new[1], (N, alpha, trial)
                assert np.array_equal(old[2], new[2]), (N, alpha, trial)


def test_sort_and_partition_agree_under_heavy_ties():
    """The tie-breaking rule is where a partition shortcut would silently diverge.

    Energies are rounded onto a tiny grid so most values are repeated, and one case is a
    single constant so EVERY element ties.  If the tail were selected by ``e <= q`` or by
    an unstable partition order, the mask -- and therefore every gradient weight -- would
    differ from the shipped code here and nowhere else.
    """
    rng = np.random.default_rng(6)
    cases = [np.round(rng.standard_normal(500), 1),        # above CVAR_SORT_CUTOFF
             np.round(rng.standard_normal(5000), 1),
             rng.integers(0, 3, 400).astype(float),        # below it
             rng.integers(0, 3, 4000).astype(float),
             np.zeros(64),
             np.zeros(2000),                               # every element ties
             np.concatenate([np.zeros(200), np.ones(200)]),
             np.concatenate([np.ones(1500), np.zeros(1500)])]
    for e in cases:
        for alpha in (0.9, 0.5, 0.25, 0.1, 0.02):
            old = qansatz.cvar(e, alpha)
            new = Q.cvar(e, alpha)
            assert old[0] == new[0] and old[1] == new[1]
            assert np.array_equal(old[2], new[2]), (alpha, e[:8])
            assert int(new[2].sum()) == max(1, int(math.ceil(alpha * e.size)))


def test_partition_is_actually_used_and_is_cheaper_where_it_is_used():
    """A correctness-preserving change that did not change the cost is not the fix.

    Above `CVAR_SORT_CUTOFF` the selection has to be materially cheaper than the sort it
    replaced, or the whole exercise was decoration.  The threshold is loose (5x against a
    measured 12-25x) because this box runs three sibling jobs and a tight bound would be a
    flake generator, but it still fails outright if the partition path secretly sorts.
    """
    import time
    rng = np.random.default_rng(8)

    def t(fn, e, reps):
        fn(e, 0.1)
        t0 = time.perf_counter()
        for _ in range(reps):
            fn(e, 0.1)
        return (time.perf_counter() - t0) / reps

    big = rng.standard_normal(1 << 16)
    old = t(qansatz.cvar, big, 20)
    new = t(Q.cvar, big, 20)
    assert new < old / 5.0, (old, new)
    # ...and below the cutoff the sort branch is taken, so it must not be SLOWER
    small = rng.standard_normal(256)
    assert t(Q.cvar, small, 2000) < t(qansatz.cvar, small, 2000) * 1.5


# ================================================== ansatz equivalence: OneLayer
def test_onelayer_ansatz_is_bit_identical_to_shipped():
    """Sampling, log p and grad log p, all ``==``. This family was not rewritten."""
    rng = np.random.default_rng(0)
    for nq in (6, 24, 64):
        for ring in (True, False):
            a0 = qansatz.OneLayerAnsatz(nq, ring=ring)
            a1 = Q.OneLayerAnsatz(nq, ring=ring)
            th = rng.uniform(0.2, math.pi - 0.2, nq)
            b0 = a0.sample(th, 300, np.random.default_rng(31))
            b1 = a1.sample(th, 300, np.random.default_rng(31))
            assert np.array_equal(b0, b1)
            assert np.array_equal(a0.logp(th, b0), a1.logp(th, b1))
            assert np.array_equal(a0.grad_logp(th, b0), a1.grad_logp(th, b1))


def test_onelayer_entropy_closed_form_matches_explicit_distribution():
    """The CNOT chain is a bijection, so the output entropy is the Bernoulli layer's.

    Checked against the explicit ``2**n`` distribution and its analytic gradient against
    finite differences, because a closed form that is merely plausible is not a closed form.
    """
    rng = np.random.default_rng(12)
    nq = 10
    an = Q.OneLayerAnsatz(nq, ring=True)
    th = rng.uniform(0.3, math.pi - 0.3, nq)
    H, g = Q.entropy_onelayer(an, th)
    p = an.probs(th)
    H_exp = float(-(p * np.log(np.maximum(p, 1e-300))).sum())
    assert abs(H - H_exp) < 1e-10, (H, H_exp)
    fd = np.zeros(nq)
    for q in range(nq):
        tp, tm = th.copy(), th.copy()
        tp[q] += 1e-6
        tm[q] -= 1e-6
        fd[q] = (Q.entropy_onelayer(an, tp)[0] - Q.entropy_onelayer(an, tm)[0]) / 2e-6
    assert round(_cos(g, fd), 6) == 1.000000, _cos(g, fd)


# ======================================================= ansatz equivalence: MPS
def _contract_padded(A):
    v = np.zeros((1, A.shape[1]))
    v[0, 0] = 1.0
    for q in range(A.shape[0]):
        v = np.einsum("db,bsr->dsr", v, A[q]).reshape(-1, A.shape[3])
    return v[:, 0]


def _contract_ragged(A):
    v = np.ones((1, 1))
    for t in A:
        v = np.einsum("db,bsr->dsr", v, t).reshape(-1, t.shape[2])
    return v[:, 0]


def test_mps_padded_build_contracts_to_the_shipped_state():
    """The padded ``(n, chi, 2, chi)`` array is the same MPS, tensor for tensor.

    The rewrite reorganised the CNOT chain from a per-site Python loop into one MPO applied
    to every site at once, which is only legitimate if the bond-index conventions survive.
    Contracting both representations to a full statevector is the check that cannot be
    fooled by a consistent relabelling that is nevertheless wrong at the boundaries -- and
    the boundaries (site 0 is never a target, site n-1 is never a control) are exactly
    where a batched MPO gets it wrong.
    """
    rng = np.random.default_rng(0)
    for n, L, fry in ((4, 1, True), (6, 2, True), (5, 1, False), (7, 2, False),
                      (6, 3, True)):
        old = qansatz.MPSAnsatz(n, L, final_ry=fry)
        new = Q.MPSAnsatz(n, L, final_ry=fry)
        th = rng.normal(0, 1.1, old.n_params())
        psi_old = _contract_ragged(old._build_np(th))
        psi_new = _contract_padded(new.build(th))
        assert np.abs(psi_old - psi_new).max() == 0.0, (n, L, fry)
        # and both are the independent dense circuit
        psi_ref = _dense_state(th, n, L + int(fry), L, ring=False)
        assert np.abs(psi_new - psi_ref).max() < 1e-12, (n, L, fry)


def test_mps_logp_matches_legacy_within_the_stated_tolerance():
    """Log p after the log-depth contraction, against the shipped left-to-right loop.

    This is the ONE place the reduction order changed, so this is the one place a tolerance
    is honest rather than lazy.  Matrix multiplication is associative, so the two compute
    the same product; they differ only in floating-point accumulation order.  The measured
    difference is quoted in the assertion so a regression that made it worse would show.
    """
    import torch
    rng = np.random.default_rng(0)
    worst = 0.0
    for n, L, fry in ((6, 1, True), (9, 2, True), (12, 1, False), (10, 2, False),
                      (7, 3, True)):
        old = qansatz.MPSAnsatz(n, L, final_ry=fry)
        new = Q.MPSAnsatz(n, L, final_ry=fry)
        th = rng.normal(0, 1.1, old.n_params())
        allb = Q.all_bitstrings(n)
        lp0 = old.logp_torch(torch.as_tensor(th),
                             torch.as_tensor(allb.astype(np.int64))).numpy()
        lp1 = new.logp(th, allb)
        worst = max(worst, float(np.abs(lp0 - lp1).max()))
        assert abs(np.exp(lp1).sum() - 1.0) < 1e-12, (n, L, fry)
    assert worst < 1e-11, worst


def test_mps_grad_logp_matches_legacy_within_the_stated_tolerance():
    rng = np.random.default_rng(1)
    worst = 0.0
    for n, L, fry in ((6, 1, True), (9, 2, True), (16, 1, True), (12, 2, False)):
        old = qansatz.MPSAnsatz(n, L, final_ry=fry)
        new = Q.MPSAnsatz(n, L, final_ry=fry)
        th = rng.normal(0, 1.1, old.n_params())
        bits = old.sample(th, 96, np.random.default_rng(3))
        w = rng.standard_normal(96)
        g0 = old.grad_logp(th, bits, w)
        g1 = new.grad_logp(th, bits, w)
        worst = max(worst, float(np.abs(g0 - g1).max()
                                 / max(np.abs(g0).max(), 1e-300)))
    assert worst < 1e-11, worst


def test_mps_sampling_is_bit_identical_to_shipped():
    """Same rng, same bits.

    The sampler is autoregressive, so a difference of one part in 1e15 in a conditional
    probability could in principle flip a bit.  It does not: the environments are contracted
    from the same tensors, and the per-site quadratic form is the same arithmetic with
    ``A`` and ``R`` folded together beforehand.
    """
    rng = np.random.default_rng(0)
    for n, L, fry, ent in ((12, 1, True, "cnot"), (20, 2, True, "cnot"),
                           (30, 1, False, "cnot"), (16, 1, False, "none")):
        old = qansatz.MPSAnsatz(n, L, final_ry=fry, entangler=ent)
        new = Q.MPSAnsatz(n, L, final_ry=fry, entangler=ent)
        th = rng.normal(0, 0.7, old.n_params())
        s0 = old.sample(th, 512, np.random.default_rng(17))
        s1 = new.sample(th, 512, np.random.default_rng(17))
        assert np.array_equal(s0, s1), (n, L, fry, ent, int((s0 != s1).sum()))


def test_mps_sampler_reproduces_the_exact_distribution():
    """The sampler draws from the state it claims to, not from a convenient proxy.

    An exact contraction that is sampled wrongly would pass every algebraic test above and
    still make the search meaningless, so the empirical histogram over the full register is
    compared against the analytic probabilities.
    """
    n = 8
    an = Q.MPSAnsatz(n, 2, final_ry=True)
    th = np.random.default_rng(2).normal(0, 1.0, an.n_params())
    p = an.probs(th)
    bits = an.sample(th, 400_000, np.random.default_rng(19))
    idx = bits.astype(np.int64) @ (1 << np.arange(n - 1, -1, -1))
    emp = np.bincount(idx, minlength=1 << n) / 400_000.0
    assert np.abs(emp - p).max() < 5e-3, np.abs(emp - p).max()
    assert np.abs(emp - p).sum() < 0.05


def test_mps_bond_dimension_is_exact_not_truncated():
    """chi == 2**layers with no SVD anywhere: the simulation is exact by construction.

    If a truncation were ever introduced to buy speed, the state would stop being the
    circuit's and every downstream number would be a plausible-looking fiction.  The bond
    dimension is asserted against the topological bound, and the state is asserted
    normalised -- a truncated MPS is not.
    """
    for n, L, fry in ((10, 1, True), (10, 2, True), (12, 3, False)):
        a = Q.MPSAnsatz(n, L, final_ry=fry)
        assert a.chi == 2 ** L, (n, L, a.chi)
        A = a.build(np.random.default_rng(0).normal(0, 1.0, a.n_params()))
        assert A.shape == (n, a.chi, 2, a.chi)
        assert abs(float(a.norm(A)) - 1.0) < 1e-12


# ============================================= statevector circuit equivalence
def test_statevector_probs_bit_identical_to_legacy():
    """`StatevectorCircuit` against `s8/integrate.Circuit`.

    ``probs`` is asserted with ``==``: the composed entangler permutation is exact and the
    RY update is the same two elementwise expressions, so there is no reason for a single
    bit to move.  The parameter-shift gradient gets a tolerance only because its final
    contraction is now one matmul instead of P separate dot products.
    """
    m8 = _s8()
    rng = np.random.default_rng(0)
    for n, L in ((5, 2), (7, 3), (10, 3)):
        c0, c1 = m8.Circuit(n, L), Q.StatevectorCircuit(n, L)
        assert c0.n_params() == c1.n_params()
        E = rng.standard_normal(c1.dim)
        for _ in range(3):
            th = rng.normal(0, 1.0, c1.n_params())
            assert np.array_equal(c0.probs(th), c1.probs(th)), (n, L)
            assert np.array_equal(c0.state(th), c1.state(th)), (n, L)
            for alpha in (0.1, 0.25):
                g0 = m8.grad_cvar_paramshift(c0, th, E, alpha)
                g1 = Q.grad_cvar_paramshift(c1, th, E, alpha)
                assert np.abs(g0 - g1).max() < 1e-14, (n, L, alpha)
                v0 = m8.cvar_exact(E, c0.probs(th), alpha)
                v1 = Q.cvar_exact(E, c1.probs(th), alpha)
                assert v0[0] == v1[0] and v0[1] == v1[1]
                assert np.array_equal(v0[2], v1[2])


def test_statevector_matches_independent_dense_simulation():
    """And both agree with the Kronecker-product reference built in this file."""
    rng = np.random.default_rng(3)
    for n, L in ((4, 1), (6, 2), (7, 3)):
        c = Q.StatevectorCircuit(n, L, ring=True)
        th = rng.normal(0, 1.0, c.n_params())
        ref = _dense_state(th, n, L, L, ring=True)
        assert np.abs(c.state(th) - ref).max() < 1e-12, (n, L)


def test_batched_probs_equals_looping_one_at_a_time():
    """Batching over parameter vectors must not couple them.

    The AMPLITUDES are bit-identical (``==``) because every gate is an elementwise
    expression on an independent row.  The normalised probabilities differ by up to 4e-17,
    and only there: ``p.sum(axis=1)`` on a (B, dim) array blocks its pairwise summation
    differently from the same reduction on one row.  That is quantified rather than hidden,
    and it never reaches the search -- the search calls `probs` with a single parameter
    vector, which `test_statevector_probs_bit_identical_to_legacy` pins at ``==``.
    """
    rng = np.random.default_rng(4)
    c = Q.StatevectorCircuit(8, 3)
    TH = rng.normal(0, 1.0, (17, c.n_params()))
    batched_s = c.states_batch(TH)
    batched_p = c.probs_batch(TH)
    for i in range(TH.shape[0]):
        assert np.array_equal(batched_s[i], c.state(TH[i])), i
        assert np.abs(batched_p[i] - c.probs(TH[i])).max() < 1e-16, i


def test_free_energy_and_run_cvar_vqe_match_legacy():
    """The free-energy objective and the driver that optimises it, against s8.

    Iterations, restarts, learning rate and seed handling are unchanged, so the optimisation
    TRAJECTORY has to land in the same place; only the per-gradient simulation was batched.
    """
    m8 = _s8()
    rng = np.random.default_rng(5)
    n, L = 7, 3
    c0, c1 = m8.Circuit(n, L), Q.StatevectorCircuit(n, L)
    E = rng.standard_normal(c1.dim)
    th = rng.normal(0, 1.0, c1.n_params())
    for T in (0.0, 0.3):
        f0, g0, p0, v0, H0 = m8.free_energy(c0, th, E, 0.25, T)
        f1, g1, p1, v1, H1 = Q.free_energy(c1, th, E, 0.25, T)
        assert f0 == f1 and v0 == v1 and H0 == H1
        assert np.array_equal(p0, p1)
        assert np.abs(g0 - g1).max() < 1e-13, np.abs(g0 - g1).max()
    for T in (0.1, 0.3):
        r0 = m8.run_cvar_vqe(E, 0.25, T, n=n, layers=L, iters=30, seed=0)
        r1 = Q.run_cvar_vqe(E, 0.25, T, n=n, layers=L, iters=30, seed=0)
        assert np.abs(r0[0] - r1[0]).max() < 1e-10, np.abs(r0[0] - r1[0]).max()
        assert abs(r0[1] - r1[1]) < 1e-10 and abs(r0[2] - r1[2]) < 1e-10


# ============================================================== PennyLane path
def test_pennylane_bound_tape_matches_the_qnode_exactly():
    """The rebuilt-tape fast path is the same PennyLane circuit on the same device.

    A QNode reconstructs its tape and re-runs the transform program on every call; this path
    builds the tape once and rebinds the angles.  Same gates, same wires, same
    ``lightning.qubit``, so the probabilities are asserted with ``==``.
    """
    import pennylane as qml
    rng = np.random.default_rng(0)
    for n, L, ring in ((4, 2, True), (9, 3, True), (12, 4, False)):
        fast = Q.build_global_circuit(n, L, ring=ring)
        dev = qml.device("lightning.qubit", wires=n)

        @qml.qnode(dev)
        def ref(params, n=n, L=L, ring=ring):
            p = np.reshape(np.asarray(params, float), (L, n))
            for l in range(L):
                for q in range(n):
                    qml.RY(float(p[l][q]), wires=q)
                for q in range(n - 1):
                    qml.CNOT(wires=[q, q + 1])
                if ring and n > 2:
                    qml.CNOT(wires=[n - 1, 0])
            return qml.probs(wires=range(n))

        for _ in range(3):
            th = rng.normal(np.pi / 2, 1.0, n * L)
            assert np.array_equal(np.asarray(fast(th)), np.asarray(ref(th))), (n, L)


def test_pennylane_circuit_matches_the_dense_reference():
    """...and PennyLane's answer is the circuit this module says it is."""
    rng = np.random.default_rng(1)
    n, L = 6, 2
    fast = Q.build_global_circuit(n, L, ring=True)
    th = rng.normal(np.pi / 2, 1.0, n * L)
    assert np.abs(np.asarray(fast(th)) - _dense_state(th, n, L, L, ring=True) ** 2
                  ).max() < 1e-12


def test_global_vqe_cvar_helpers_match_shipped():
    """``cvar_from_samples`` is the objective SPSA sees, so it is asserted BIT-identical.

    One ulp is not a rounding detail here: SPSA takes a finite difference of two noisy
    evaluations and divides by a small ``ck``, so a last-bit change in the objective can
    flip the sign of a step and fork the whole trajectory.  That is why this routine sorts
    its selected values before averaging instead of taking the mean over a boolean mask --
    the shipped code summed in ascending order and so does this one.
    """
    rng = np.random.default_rng(2)
    for N in (16, 257, 2048, 8192):
        for e in (rng.standard_normal(N), np.round(rng.standard_normal(N), 1),
                  np.zeros(N)):
            for alpha in (0.05, 0.15, 0.5, 1.0):
                assert VQE.cvar_from_samples(e, alpha) == Q.cvar_from_samples(e, alpha), \
                    (N, alpha)
        e = rng.standard_normal(N)
        p = rng.random(N)
        for alpha in (0.05, 0.15, 1.0):
            assert abs(VQE.cvar_from_distribution(e, p, alpha)
                       - Q.cvar_from_distribution(e, p, alpha)) < 1e-12
    assert Q.n_parameters(11, 4) == VQE.n_parameters(11, 4)


# ================================================ 4. old vs new: the search itself
class _ToyRep:
    """A torsion representation with the interface the driver uses and nothing else.

    Deterministic and self-contained: `run_vqe` needs ``bits_per_residue``, ``n_residues``,
    ``n_states`` and the ``_phi`` / ``_psi`` tables, and using a real target here would make
    the trajectory comparison depend on the distogram, the ESM bank and the peptide
    database -- none of which this module owns or is testing.
    """

    def __init__(self, n_res=9, bits=2, seed=0):
        self.n_residues = int(n_res)
        self.bits_per_residue = int(bits)
        self.n_states = 1 << int(bits)
        rng = np.random.default_rng(seed)
        self._phi = rng.uniform(-np.pi, np.pi, (self.n_residues, self.n_states))
        self._psi = rng.uniform(-np.pi, np.pi, (self.n_residues, self.n_states))


class _ToyObjective:
    """A deterministic Potts-like score plus a CA trace. No randomness at call time."""

    def __init__(self, rep, seed=0):
        self.rep = rep
        self.n = rep.n_residues
        self.mrf = None
        self.w_mrf = 0.0
        self.w_clash = 0.0
        rng = np.random.default_rng(seed)
        k = rep.n_states
        self.h = rng.standard_normal((self.n, k))
        self.J = rng.standard_normal((self.n - 1, k, k))
        d = rng.standard_normal((k, 3))
        self.dirs = d / np.linalg.norm(d, axis=1, keepdims=True)
        self.n_calls = 0
        self.n_structures = 0

    def __call__(self, states):
        S = np.asarray(states, int)
        e = self.h[np.arange(self.n)[None, :], S].sum(1)
        e += self.J[np.arange(self.n - 1)[None, :], S[:, :-1], S[:, 1:]].sum(1)
        ca = np.cumsum(self.dirs[S] * 3.8, axis=1)
        self.n_calls += 1
        self.n_structures += len(S)
        return e, ca


def _run_pair(kind, iters=25, shots=192, seed=0, update="cvar_grad"):
    """One restart through the shipped driver and through `core.quantum`, same seed."""
    rep = _ToyRep(seed=1)
    nq = rep.bits_per_residue * rep.n_residues
    o0, o1 = _ToyObjective(rep, seed=2), _ToyObjective(rep, seed=2)
    if kind == "chain":
        a0 = qansatz.OneLayerAnsatz(nq, ring=True)
        a1 = Q.OneLayerAnsatz(nq, ring=True)
    else:
        a0 = qansatz.MPSAnsatz(nq, 1, final_ry=True)
        a1 = Q.MPSAnsatz(nq, 1, final_ry=True)
    th0 = np.full(a0.n_params(), np.pi / 2) + \
        np.random.default_rng(99).normal(0, 0.35, a0.n_params())
    tr0, tr1 = [], []
    r0 = FV.run_vqe(o0, a0, th0, iters, shots, np.random.default_rng([seed, 0]),
                    update=update, trace=tr0)
    r1 = Q.run_vqe(o1, a1, th0, iters, shots, np.random.default_rng([seed, 0]),
                   update=update, baseline="tail", trace=tr1)
    return r0, tr0, r1, tr1


def test_search_trajectory_matches_shipped_driver_on_the_chain_ansatz():
    """Old vs new, the whole optimisation trajectory, bit-for-bit.

    `foldvqe.run_vqe` and `core.quantum.run_vqe` are handed the same objective, the same
    initial angles and the same rng, with the new driver forced onto the SHIPPED baseline so
    the only thing being compared is the rewrite.  On the one-layer chain nothing about the
    arithmetic changed, so every CVaR value, every best-so-far and the final angle vector
    are asserted with ``==``.
    """
    for update in ("cvar_grad", "moment"):
        r0, tr0, r1, tr1 = _run_pair("chain", update=update)
        assert len(tr0) == len(tr1) == 25
        for a, b in zip(tr0, tr1):
            assert a["cvar"] == b["cvar"] and a["min"] == b["min"], (update, a, b)
            assert a["best"] == b["best"] and a["alpha"] == b["alpha"]
        assert np.array_equal(r0["theta"], r1["theta"]), update
        assert r0["best_energy"] == r1["best_energy"], update
        assert np.array_equal(r0["best_states"], r1["best_states"]), update


def test_search_trajectory_matches_shipped_driver_on_the_mps_ansatz():
    """Same comparison on the MPS family, where the reduction order DID change.

    The log-depth contraction moves log p by ~1e-13, which the optimiser then amplifies
    across 25 Adam steps.  That is quantified rather than asserted away: the trajectory is
    required to agree to 1e-8 in the objective and the final angles to 1e-8, which is
    several orders tighter than any quantity this search reports.
    """
    r0, tr0, r1, tr1 = _run_pair("mps")
    dc = max(abs(a["cvar"] - b["cvar"]) for a, b in zip(tr0, tr1))
    db = max(abs(a["best"] - b["best"]) for a, b in zip(tr0, tr1))
    dth = float(np.abs(r0["theta"] - r1["theta"]).max())
    assert dc < 1e-8, dc
    assert db < 1e-8, db
    assert dth < 1e-8, dth


def test_new_default_baseline_changes_the_search_and_that_is_the_point():
    """The corrected gradient must actually alter the trajectory.

    If ``baseline="const"`` and ``baseline="tail"`` produced the same search, the fix would
    be cosmetic and the audit above would be measuring nothing that matters.  They do not.
    """
    rep = _ToyRep(seed=1)
    nq = rep.bits_per_residue * rep.n_residues
    a = Q.OneLayerAnsatz(nq, ring=True)
    th0 = np.full(nq, np.pi / 2)
    outs = {}
    for b in ("const", "tail"):
        outs[b] = Q.run_vqe(_ToyObjective(rep, seed=2), a, th0, 25, 192,
                            np.random.default_rng([0, 0]), baseline=b)
    assert not np.allclose(outs["const"]["theta"], outs["tail"]["theta"], atol=1e-6)


def test_reservoir_replay_is_identical_to_sequential_offering():
    """Concurrent restarts: the record-and-replay contract.

    The Reservoir is order-dependent -- an offer inside the diversity radius overwrites its
    neighbour -- so restarts cannot share one.  They do not need to: a restart's offers
    depend only on its own sampling.  This runs the restarts sequentially into a shared
    reservoir, then again into per-restart logs replayed in restart order, and asserts the
    two reservoirs are the same object down to the last float.
    """
    rep = _ToyRep(seed=1)
    nq = rep.bits_per_residue * rep.n_residues
    a = Q.OneLayerAnsatz(nq, ring=True)
    th0 = np.full(nq, np.pi / 2)
    o = _ToyObjective(rep, seed=2)

    direct = Q.Reservoir(radius=1.2, capacity=16)
    for r in range(4):
        Q.run_vqe(o, a, th0, 12, 128, np.random.default_rng([0, r]),
                  reservoir=direct, random_frac=0.05)

    logs = [Q._OfferLog() for _ in range(4)]
    for r in range(4):
        Q.run_vqe(o, a, th0, 12, 128, np.random.default_rng([0, r]),
                  reservoir=logs[r], random_frac=0.05)
    replayed = Q.Reservoir(radius=1.2, capacity=16)
    for lg in logs:
        lg.replay(replayed)

    assert direct.scores == replayed.scores
    assert len(direct.cas) == len(replayed.cas) > 1
    for x, y in zip(direct.cas, replayed.cas):
        assert np.array_equal(x, y)
    for x, y in zip(direct.states, replayed.states):
        assert np.array_equal(x, y)


def test_reservoir_and_selection_match_shipped():
    """`Reservoir` and `basin_select` are ports, so they get equality assertions."""
    rng = np.random.default_rng(0)
    r0, r1 = FV.Reservoir(radius=1.5, capacity=8), Q.Reservoir(radius=1.5, capacity=8)
    for _ in range(12):
        s = rng.standard_normal(6)
        ca = rng.standard_normal((6, 10, 3)) * 4.0
        st = rng.integers(0, 4, (6, 10))
        r0.offer_batch(s, ca, st)
        r1.offer_batch(s, ca, st)
    assert r0.scores == r1.scores
    for a, b in zip(r0.cas, r1.cas):
        assert np.array_equal(a, b)
    p0, i0 = FV.basin_select(r0.scores, r0.cas, radius=2.0)
    p1, i1 = Q.basin_select(r1.scores, r1.cas, radius=2.0)
    assert p0 == p1 and i0["members"] == i1["members"]
    assert i0["basin_size"] == i1["basin_size"]


def test_warm_start_and_marginal_inversion_match_shipped():
    """The closed-form marginal inversion is the warm start; it is a port, so ``==``.

    It is also checked to BE an inversion, on marginals the plain chain can actually reach.
    That qualifier is the documented constraint, not a convenience: ``b_q`` is the prefix
    XOR of independent Bernoulli bits, so ``|1 - 2 P(b_q = 1)|`` is a running product and
    must be non-increasing along the register.  Reachable marginals are therefore generated
    FROM a Bernoulli vector rather than drawn at random, and the second half of this test
    shows that random marginals are genuinely unreachable -- which is the measurement that
    makes the trailing RY layer of `MPSAnsatz` load-bearing rather than decorative.
    """
    rng = np.random.default_rng(0)
    for nq in (6, 12, 20):
        s = rng.uniform(0.1, 0.9, nq)                      # a reachable configuration
        m = (1.0 - np.cumprod(1.0 - 2.0 * s)) / 2.0
        assert np.array_equal(FV.marginals_to_angles(m), Q.marginals_to_angles(m))
        th = Q.marginals_to_angles(m)
        # The inversion is asserted on the MARGINALS, which is what it promises and what
        # the warm start uses. It does NOT recover the underlying Bernoulli parameters to
        # the same accuracy at long registers, and that is a property of the target rather
        # than of the arithmetic: |1 - 2 P(b_q = 1)| is a running product, so by q=30 it has
        # decayed to 2e-14 and every remaining marginal is 0.5 to within round-off. Several
        # different s vectors then produce the same marginals.
        back = (1.0 - np.cumprod(1.0 - 2.0 * np.sin(th / 2.0) ** 2)) / 2.0
        assert np.abs(back - m).max() < 1e-12, (nq, np.abs(back - m).max())
        if nq <= 12:
            assert np.abs(np.sin(th / 2.0) ** 2 - s).max() < 1e-9, nq
        b = Q.OneLayerAnsatz(nq, ring=False).sample(th, 200_000,
                                                    np.random.default_rng(1))
        assert np.abs(b.mean(0) - m).max() < 0.01, np.abs(b.mean(0) - m).max()

    # the constraint, measured: an arbitrary marginal vector is NOT reachable by the chain
    bad = rng.uniform(0.05, 0.95, 8)
    got = np.sin(Q.marginals_to_angles(bad) / 2.0) ** 2
    reached = (1.0 - np.cumprod(1.0 - 2.0 * got)) / 2.0
    assert np.abs(reached - bad).max() > 0.1, np.abs(reached - bad).max()

    pi = rng.random((7, 8))
    pi /= pi.sum(1, keepdims=True)
    assert np.array_equal(FV.state_marginals_to_bits(pi, 3),
                          Q.state_marginals_to_bits(pi, 3))


def test_marginal_inversion_degenerates_identically_to_shipped_at_p_equals_half():
    """A shipped edge case, ported rather than silently repaired.

    When a marginal is exactly 0.5 the running product ``1 - 2 P(b_q = 1)`` is zero and the
    next ratio is 0/0, so the inversion returns nan from that wire on.  It happens for a
    long chain of interior Bernoulli parameters, where the product underflows to exactly
    zero.  Real prior marginals do not sit on 0.5 to machine precision, so this has never
    fired in production -- but the consolidated copy must reproduce it rather than quietly
    guard it, because a guard would change the warm start on inputs the shipped code
    rejects.  Both implementations produce nan in exactly the same positions.
    """
    s = np.random.default_rng(0).uniform(0.1, 0.9, 60)
    m = (1.0 - np.cumprod(1.0 - 2.0 * s)) / 2.0
    with np.errstate(invalid="ignore"):
        a, b = FV.marginals_to_angles(m), Q.marginals_to_angles(m)
    assert np.array_equal(np.isnan(a), np.isnan(b))
    assert np.isnan(a).any()
    assert np.array_equal(a[~np.isnan(a)], b[~np.isnan(b)])


# ================================================================== genuineness
def test_the_search_does_not_enumerate():
    """The optimiser must touch a vanishing fraction of the register.

    "Genuine VQE" is falsifiable and this is the falsification: at 18 qubits the state space
    is 262,144, and a search that quietly enumerated it would show up here as a bitstring
    count of that order.  A real sampled VQE at 25 iterations x 192 shots can touch at most
    4,800 -- under 2% -- and in practice far fewer once the distribution concentrates.
    """
    rep = _ToyRep(n_res=9, bits=2, seed=1)
    nq = rep.bits_per_residue * rep.n_residues
    assert nq == 18
    seen = set()
    an = Q.OneLayerAnsatz(nq, ring=True)

    class Counting(_ToyObjective):
        def __call__(self, states):
            for s in np.asarray(states, int):
                seen.add(s.tobytes())
            return super().__call__(states)

    Q.run_vqe(Counting(rep, seed=2), an, np.full(nq, np.pi / 2), 25, 192,
              np.random.default_rng([0, 0]))
    assert len(seen) <= 25 * 192
    assert len(seen) < 0.02 * (1 << nq), (len(seen), 1 << nq)


def test_enumeration_helpers_refuse_to_be_load_bearing():
    """The verification instruments are capped so they cannot drift into the search path."""
    import pytest
    with pytest.raises(ValueError):
        Q.all_bitstrings(21)
    assert Q.all_bitstrings(4).shape == (16, 4)


def test_cvar_rejects_degenerate_alpha_and_empty_samples():
    import pytest
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            Q.cvar(np.arange(10.0), bad)
    with pytest.raises(ValueError):
        Q.cvar(np.array([]), 0.5)
    with pytest.raises(ValueError):
        Q.cvar_gradient(Q.OneLayerAnsatz(4), np.zeros(4),
                        np.zeros((2, 4), np.uint8), np.zeros(2), 0.5,
                        baseline="mean-of-tail")


def test_objective_and_hamiltonian_are_the_shipped_classes():
    """`FoldObjective` and `FoldingHamiltonian` are ports; their signatures must not drift.

    They define the energy the VQE minimises, so a silent change to either would invalidate
    every comparison the project has made.  Checked structurally (the same public methods
    and the same defaults) because instantiating them needs the distogram and the peptide
    database, which belong to sibling modules rather than to this one.
    """
    import inspect
    import hamiltonian as HAM
    import objective as OBJ
    for new, old in ((Q.FoldObjective, OBJ.FoldObjective),
                     (Q.FoldingHamiltonian, HAM.FoldingHamiltonian)):
        s_new = inspect.signature(new.__init__)
        s_old = inspect.signature(old.__init__)
        assert list(s_new.parameters) == list(s_old.parameters), new
        for k, p in s_old.parameters.items():
            assert s_new.parameters[k].default == p.default, (new, k)
        for name in dir(old):
            if name.startswith("_") or not callable(getattr(old, name, None)):
                continue
            assert hasattr(new, name), (new, name)
    assert Q.MAX_NATIVE_CHI_SCAN_BITS == HAM.MAX_NATIVE_CHI_SCAN_BITS
    assert Q.MIN_USEFUL_SPSA_ITERS == VQE.MIN_USEFUL_SPSA_ITERS


def test_geometry_backend_agrees_with_the_root_module():
    """`FoldObjective` builds through `core.backend("geometry")`, so that has to agree.

    The switch is the project's design and this module uses it rather than importing
    `protein_geometry` directly -- which means a consolidated `core.geometry` going live
    silently changes what `fold` scores.  This asserts, on exactly the call `FoldObjective`
    makes, that whichever backend is live returns the same coordinates as the root module.
    It is a guard on someone else's file, but the consequence lands here.
    """
    import protein_geometry as root
    if Q.geo is root:
        return                                       # legacy mode: nothing to compare
    rng = np.random.default_rng(0)
    phi = rng.uniform(-np.pi, np.pi, (32, 14))
    psi = rng.uniform(-np.pi, np.pi, (32, 14))
    a, b = Q.geo.build_backbone_batch(phi, psi), root.build_backbone_batch(phi, psi)
    assert set(a) == set(b)
    for k in a:
        assert np.abs(np.asarray(a[k]) - np.asarray(b[k])).max() < 1e-10, k
    ref = rng.standard_normal((14, 3)) * 4.0
    assert np.abs(Q.geo.ca_rmsd_batch(np.asarray(a["CA"]), ref)
                  - root.ca_rmsd_batch(np.asarray(b["CA"]), ref)).max() < 1e-10


def test_spsa_is_the_shipped_optimiser():
    """SPSA's arithmetic is untouched: same trajectory on a deterministic objective."""
    def f(x):
        return float(np.sum((x - 0.3) ** 2) + 0.1 * np.sum(np.sin(3 * x)))

    x0 = np.linspace(-1, 1, 11)
    r0 = VQE._spsa(f, x0, 60, np.random.default_rng(7))
    r1 = Q._spsa(f, x0, 60, np.random.default_rng(7))
    assert np.array_equal(r0.x, r1.x)
    assert r0.fun == r1.fun


# ================================================================ fallback runner
def _main():
    fns = [(k, v) for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    fails = []
    for name, fn in fns:
        try:
            fn()
            print(f"PASS  {name}")
        except Exception as exc:
            fails.append(name)
            print(f"FAIL  {name}   {type(exc).__name__}: {exc}")
    print(f"\n{len(fns) - len(fails)}/{len(fns)} passed"
          + (f"   FAILURES: {fails}" if fails else ""))
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(_main())
