"""s26/q_tests.py -- unit tests for s26/q_adapt.py, s26/q_dla.py, s26/q_var.py on SYNTHETIC
energies at n <= 6.  No native, no RMSD, no pipeline.  Runs under pytest or as a script:

    python s26/q_tests.py
    python -m pytest s26/q_tests.py -q

The checks the brief names: the simulator against core.quantum.StatevectorCircuit for the
fixed ansatz and against pennylane default.qubit for Pauli rotations, both to 1e-12; the
parameter-shift gradients against finite differences; the pool-gradient lemma; the symbolic DLA
closure against the numerical route; and the ADAPT loop's bookkeeping.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s26 import q_adapt as QA                                             # noqa: E402
from s26 import q_dla as QD                                               # noqa: E402


def _rand_E(dim, rng):
    E = rng.normal(0, 1, dim)
    return (E - E.mean()) / E.std()


# ------------------------------------------------------------- Pauli action
def test_pauli_action_matches_dense_matrices():
    """(-iP) psi from the bit tricks equals the Kronecker-product matrix, n = 2..4."""
    rng = np.random.default_rng(0)
    for n in (2, 3, 4):
        for _ in range(30):
            word = "".join(rng.choice(list("IXYZ"), n))
            P = QA.Pauli.from_word(word)
            M = (-1j) * QD.dense_pauli(P)
            psi = rng.normal(size=1 << n) + 1j * rng.normal(size=1 << n)
            perm, coef = QA.pauli_action(P)
            got = coef * psi[perm]
            assert np.abs(got - M @ psi).max() < 1e-13, word
            assert np.isrealobj(coef) == (P.n_y() % 2 == 1), word


def test_pauli_word_roundtrip_and_commutation():
    rng = np.random.default_rng(1)
    for n in (3, 5):
        for _ in range(50):
            w = "".join(rng.choice(list("IXYZ"), n))
            P = QA.Pauli.from_word(w)
            assert P.word() == w
            Q = QA.Pauli.from_word("".join(rng.choice(list("IXYZ"), n)))
            A, B = QD.dense_pauli(P), QD.dense_pauli(Q)
            c = np.abs(A @ B - B @ A).max() < 1e-12
            assert P.commutes(Q) == c
            if not c:
                prod = QD.dense_pauli(P.times(Q))
                R = A @ B
                # the product is the string up to a phase
                k = np.flatnonzero(np.abs(prod.ravel()) > 0.5)[0]
                ph = R.ravel()[k] / prod.ravel()[k]
                assert np.abs(R - ph * prod).max() < 1e-12


# ------------------------------------------------------------- the simulator
def test_pauli_rotations_match_pennylane_default_qubit():
    """Products of exp(-i theta P/2) on |0..0> against pennylane, n = 2..6, to 1e-12."""
    import pennylane as qml
    rng = np.random.default_rng(2)
    for n in (2, 3, 4, 5, 6):
        for trial in range(4):
            k = int(rng.integers(3, 9))
            words = []
            for _ in range(k):
                w = "".join(rng.choice(list("IXYZ"), n))
                if w == "I" * n:
                    w = "Y" + "I" * (n - 1)
                words.append(w)
            ops = [QA.Pauli.from_word(w) for w in words]
            th = rng.normal(0, 1.2, k)
            circ = QA.PauliCircuit(n, ops)
            psi = circ.state(th)
            dev = qml.device("default.qubit", wires=n)

            @qml.qnode(dev)
            def qn():
                for w, t in zip(words, th):
                    qml.PauliRot(t, w, wires=range(n))
                return qml.state()

            ref = np.asarray(qn(), complex)
            assert np.abs(np.asarray(psi, complex) - ref).max() < 1e-12, (n, words)
            assert np.abs(circ.probs(th) - np.abs(ref) ** 2).max() < 1e-12
            if all(P.n_y() % 2 == 1 for P in ops):
                assert circ.real and np.isrealobj(psi)


def test_fixed_ansatz_rebuilt_matches_core_statevector():
    """FixedAnsatzAlt (RY as Pauli rotations + an independently composed permutation)
    against core.quantum.StatevectorCircuit, n = 3..7, L = 1..3, to 1e-12."""
    from core import quantum as Q
    rng = np.random.default_rng(3)
    for n in (3, 4, 5, 6, 7):
        for L in (1, 2, 3):
            c0 = Q.StatevectorCircuit(n, L)
            c1 = QA.FixedAnsatzAlt(n, L)
            for _ in range(3):
                th = rng.normal(0, 0.6, n * L)
                assert np.abs(c0.state(th) - c1.states_batch(th[None])[0]).max() < 1e-12
                assert np.abs(c0.probs(th) - c1.probs(th)).max() < 1e-12
    assert np.array_equal(QA.entangler_perm(7), Q.StatevectorCircuit(7, 3)._perm)


def test_batched_states_equal_single_states():
    rng = np.random.default_rng(4)
    n = 5
    ops = QA.pool_L2(n)[:9]
    circ = QA.PauliCircuit(n, ops)
    TH = rng.normal(0, 1.0, (7, len(ops)))
    S = circ.states_batch(TH)
    for i in range(7):
        assert np.array_equal(S[i], circ.state(TH[i]))


# ------------------------------------------------------------- gradients
def test_free_energy_grad_matches_finite_differences():
    rng = np.random.default_rng(5)
    n = 5
    ops = QA.ry_layer_ops(n) + QA.pool_V(n)
    circ = QA.PauliCircuit(n, ops)
    E = _rand_E(circ.dim, rng)
    for alpha, T in ((1.0, 0.3), (0.25, 0.3), (0.1, 0.0), (1.0, 0.0)):
        th = rng.normal(0, 0.8, circ.n_params())
        F, g, p, v, H = QA.free_energy_grad(circ, th, E, alpha, T)
        h = 1e-5
        fd = np.empty_like(g)
        for k in range(len(th)):
            tp, tm = th.copy(), th.copy()
            tp[k] += h
            tm[k] -= h
            fd[k] = (QA.free_energy_grad(circ, tp, E, alpha, T)[0]
                     - QA.free_energy_grad(circ, tm, E, alpha, T)[0]) / (2 * h)
        cos = g @ fd / (np.linalg.norm(g) * np.linalg.norm(fd))
        assert cos > 1 - 1e-9, (alpha, T, cos)
        assert np.linalg.norm(g - fd) / np.linalg.norm(fd) < 1e-6
        # the objective is the deployed one: same value as core.quantum.free_energy's formula
        from core import quantum as Q
        v2, q2, dp2 = Q.cvar_exact(E, p, alpha)
        lp = np.log(np.maximum(p, 1e-15))
        assert abs(F - (v2 - T * float(-(p * lp).sum()))) < 1e-12


def test_pool_gradient_lemma_and_shift_rule():
    """dF/dphi|0 by the elementwise formula == literal +-pi/2 shifted circuits == FD."""
    rng = np.random.default_rng(6)
    n = 4
    circ = QA.PauliCircuit(n, QA.ry_layer_ops(n) + QA.pool_G(n)[:2])
    E = _rand_E(circ.dim, rng)
    pool = QA.pool_L2(n)
    for alpha, T in ((1.0, 0.3), (0.25, 0.3), (0.1, 0.0)):
        th = rng.normal(0, 0.7, circ.n_params())
        g1 = QA.pool_gradients(circ, th, E, alpha, T, pool)
        g2 = QA.pool_gradients_shift(circ, th, E, alpha, T, pool)
        assert np.abs(g1 - g2).max() < 1e-12
        h = 1e-6
        for k in (0, 5, len(pool) - 1):
            c2 = QA.PauliCircuit(n, circ.ops + [pool[k]])
            fp = QA.free_energy_grad(c2, np.r_[th, h], E, alpha, T)[0]
            fm = QA.free_energy_grad(c2, np.r_[th, -h], E, alpha, T)[0]
            assert abs((fp - fm) / (2 * h) - g1[k]) < 1e-6


def test_all_pool_gradients_vanish_at_the_computational_basis_state():
    """The lemma behind starting from an RY layer: from |0..0> every gradient is exactly 0."""
    rng = np.random.default_rng(7)
    n = 5
    circ = QA.PauliCircuit(n, [])
    E = _rand_E(circ.dim, rng)
    for pool in (QA.pool_V(n), QA.pool_G(n), QA.pool_L2(n), QA.pool_L2(n, complex_pool=True)):
        g = QA.pool_gradients(circ, np.zeros(0), E, 0.25, 0.3, pool)
        assert np.abs(g).max() == 0.0


# ------------------------------------------------------------- pools
def test_pool_sizes_and_y_parity():
    for n in (2, 3, 4, 5, 7):
        V, G = QA.pool_V(n), QA.pool_G(n)
        assert len(V) == 2 * n - 2 and len(G) == 2 * n - 2
        assert len(set(V)) == len(V) and len(set(G)) == len(G)
        assert all(P.n_y() % 2 == 1 for P in V + G)
        L2 = QA.pool_L2(n)
        assert len(L2) == 2 * n * n - n and len(set(L2)) == len(L2)
        assert all(P.n_y() % 2 == 1 and P.weight() <= 2 for P in L2)
        C = QA.pool_L2(n, complex_pool=True)
        assert len(C) == 3 * n + 9 * n * (n - 1) // 2
    # Tang et al.'s printed n = 3 example: {Z3 Z2 Y1, Z3 Y2, Y3, Y2}  (qubit k -> wire k-1)
    assert sorted(P.word() for P in QA.pool_V(3)) == sorted(["YZZ", "IYZ", "IIY", "IYI"])
    assert sorted(P.word() for P in QA.pool_G(3)) == sorted(["YZI", "IYZ", "IYI", "IIY"])


# ------------------------------------------------------------- ADAPT bookkeeping
def test_run_adapt_grows_records_and_lowers_F():
    rng = np.random.default_rng(8)
    n = 4
    E = _rand_E(1 << n, rng)
    for pool in ("V", "L2"):
        for opt in ("adam_best", "lbfgs", "adam_last"):
            r = QA.run_adapt(E, 0.25, 0.3, n, pool, max_params=3 * n, adam_steps=15, seed=0,
                             record_at=(4, 8, 12), eps=1e-9, optimiser=opt)
            assert r["stopped"] == "max_params" and len(r["theta"]) == 3 * n
            assert len(r["ops"]) == 3 * n and len(r["sequence"]) == 2 * n
            assert set(r["snapshots"]) == {4, 8, 12}
            for P, s in r["snapshots"].items():
                assert abs(s["p"].sum() - 1.0) < 1e-12 and s["P"] == P and not s["stopped_early"]
            F = [t["F"] for t in r["trace"]]
            if opt != "adam_last":
                # a new angle enters at 0 (F unchanged) and the re-optimiser returns a point no
                # worse than its start, so the trace is non-increasing by construction
                assert all(b <= a + 1e-12 for a, b in zip(F, F[1:])), (pool, opt, F)
            assert r["total_adam_steps"] == 15 * (1 + 2 * n)
            assert r["n_distinct_ops"] + r["consecutive_repeats"] <= 2 * n + r["n_distinct_ops"]
    r = QA.run_adapt(E, 1.0, 0.3, n, "V", max_params=3 * n, adam_steps=200, seed=0, eps=10.0)
    assert r["stopped"] == "eps" and len(r["theta"]) == n
    assert set(r["snapshots"]) == {7, 14, 21} and all(s["stopped_early"] for s in r["snapshots"].values())


def test_gibbs_state_of_the_deployed_ladder_is_a_product_state():
    """E affine in the register index -> exp(-E/T) factorises over the bits.  Measured."""
    for n in (4, 7, 10):
        dim = 1 << n
        r = np.arange(1, dim + 1, dtype=float)
        E = (r - r.mean()) / r.std()
        for T in (0.1, 0.3, 1.0):
            d = QA.product_diagnostics(E, T, n)
            assert d["kl_gibbs_to_product"] < 1e-12 and d["kl_product_to_gibbs"] < 1e-12, (n, T, d)
    # and a generic (non-affine) E is NOT a product state
    rng = np.random.default_rng(10)
    E = _rand_E(128, rng)
    assert QA.product_diagnostics(E, 0.3, 7)["kl_gibbs_to_product"] > 1e-3
    # the trained single RY layer reaches the product Gibbs state at alpha = 1 (n = 5)
    n = 5
    r = np.arange(1, (1 << n) + 1, dtype=float)
    E = (r - r.mean()) / r.std()
    a = QA.run_adapt(E, 1.0, 0.3, n, "V", max_params=3 * n, adam_steps=50, seed=0, eps=1e-3,
                     optimiser="lbfgs", record_at=(n,))
    p = a["snapshots"][n]["p"]
    assert QA.kl(p, QA.gibbs(E, 0.3)) < 1e-6, QA.kl(p, QA.gibbs(E, 0.3))
    assert a["stopped"] == "eps"


def test_matched_entropy_random_weights():
    for H_bits in (0.5, 3.0, 5.66, 6.9):
        w = QA.random_weights_matched_entropy(H_bits * math.log(2.0), 128, seed=17)
        assert abs(QA.entropy_nats(w) / math.log(2.0) - H_bits) < 1e-6
        assert abs(w.sum() - 1.0) < 1e-12 and (w >= 0).all()
    a = QA.random_weights_matched_entropy(3.0, 128, seed=1)
    b = QA.random_weights_matched_entropy(3.0, 128, seed=1)
    c = QA.random_weights_matched_entropy(3.0, 128, seed=2)
    assert np.array_equal(a, b) and not np.allclose(a, c)


def test_match_temperature_and_distinct_states():
    rng = np.random.default_rng(11)
    E = _rand_E(128, rng)
    r = np.arange(1, 129, dtype=float)
    Ez = (r - r.mean()) / r.std()
    H = QA.entropy_nats(QA.gibbs(Ez, 0.3))
    T = QA.match_temperature(E, H)
    assert abs(QA.entropy_nats(QA.gibbs(E, T)) - H) < 1e-7
    P = [QA.gibbs(Ez, 0.3), QA.gibbs(Ez, 0.3) * (1 + 1e-4 * rng.random(128)),
         QA.gibbs(E, 0.3), QA.gibbs(Ez, 1.0)]
    k, assign = QA.distinct_states(P, 0.01)
    assert k == 3 and assign == [0, 0, 1, 2]


def test_energy_variants_are_monotone_and_zrank_is_the_deployed_E():
    from core.pipeline import _zrank
    rng = np.random.default_rng(9)
    s = np.sort(rng.normal(0, 1, 128))
    s[10:14] = s[10]                              # a tie block, as the real score has
    EV = QA.energy_variants(s)
    assert np.array_equal(EV["zrank"], _zrank(s))
    for k, v in EV.items():
        assert (np.diff(v) >= -1e-12).all(), k
        assert abs(v.mean()) < 1e-9 and abs(v.std() - 1.0) < 1e-9, k


# ------------------------------------------------------------- DLA
def test_symbolic_conjugation_matches_dense_conjugation():
    """C^k Y_q C^-k as a string, checked against dense matrices at n = 4 (up to sign)."""
    n = 4
    dim = 1 << n
    C = np.eye(dim)[QA.entangler_perm(n)]           # the gather as a matrix: psi' = psi[perm]
    for q in range(n):
        Y = QA.ry_layer_ops(n)[q]
        for k in (1, 2, 3):
            P = QD.conj_entangler(Y, power=k)
            Ck = np.linalg.matrix_power(C, k)
            lhs = Ck @ QD.dense_pauli(Y) @ Ck.T
            rhs = QD.dense_pauli(P)
            ok = np.abs(lhs - rhs).max() < 1e-12 or np.abs(lhs + rhs).max() < 1e-12
            assert ok, (q, k, P.word())
            assert P.n_y() % 2 == 1


def test_dla_symbolic_equals_numeric_and_left_right_conventions_agree():
    for n in (3, 4):
        for L in (1, 2, 3):
            gr = QD.fixed_generators(n, L, convention="right")
            gl = QD.fixed_generators(n, L, convention="left")
            sym = QD.lie_closure(gr)["dim"]
            assert QD.lie_closure(gl)["dim"] == sym
            num = QD.numeric_closure(gr)["dim"]
            assert num == sym, (n, L, sym, num)
            assert sym <= QD.dim_so(n)
    # depth 1 is the abelian algebra of the n commuting RY generators
    assert QD.lie_closure(QD.fixed_generators(5, 1))["dim"] == 5
    # a known closure: {Y1, Z1 Y2} on 2 qubits generates so(4) (dim 6)?  measured, not assumed
    two = [QA.Pauli.from_word("YI"), QA.Pauli.from_word("ZY")]
    assert QD.lie_closure(two)["dim"] == QD.numeric_closure(two)["dim"]


def test_dla_closure_cap_and_pool_completeness_at_small_n():
    for n in (2, 3, 4):
        for pn in ("V", "G"):
            c = QD.lie_closure(QA.make_pool(pn, n))
            assert c["all_odd_y"] and c["dim"] <= QD.dim_so(n)
            assert c["dim"] == QD.numeric_closure(QA.make_pool(pn, n))["dim"]
    c = QD.lie_closure(QA.pool_L2(6), cap=100)
    assert c["exceeded_cap"]


# ------------------------------------------------------------- runner
def _main():
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    fails = []
    for name, fn in fns:
        try:
            fn()
            print(f"PASS  {name}", flush=True)
        except Exception as exc:                                             # noqa: BLE001
            fails.append(name)
            import traceback
            print(f"FAIL  {name}   {type(exc).__name__}: {exc}", flush=True)
            traceback.print_exc()
    print(f"\n{len(fns) - len(fails)}/{len(fns)} passed" + (f"   FAILURES: {fails}" if fails else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(_main())
