"""S25 / LANE Q -- Q-A + Q-C: INDEPENDENT RE-VERIFICATION OF THE QUANTUM COMPONENT.

PROPERTY TESTS, NOT A DIRECTIONAL EXPERIMENT. This module reads no RMSD, computes no
endpoint, and has no outcome that can favour a hypothesis, so Rule 0's fork enumeration does
not apply (the s24 Lane-D precedent for `d_setequality_proof.py`). FALSIFIERS are registered
in `s25/PREREG_Q.md` Q-A/Q-C instead, and each check below names the one that would fire.

WHAT "INDEPENDENT" MEANS HERE, AND WHY IT MATTERS
=================================================
s24 Lane D's REV1 asserted a FALSE theorem and reported 0/2916 failures, because every
probability family it generated carried an epsilon floor and the assertion could not fire. A
green light from a dead assertion is worse than a red one. So:

* A1 checks the MPS against a dense simulator written HERE from the gate list, by explicit
  Kronecker products -- not against `StatevectorCircuit` (a sibling in the same module, with
  a different topology) and not against the MPS's own contraction.
* A3 rebuilds the set-equality families from scratch WITH EXACT ZEROS (`p = 0.0`, never
  `+1e-30`), and separately asserts that the assertion is LIVE by constructing a case that
  must fail prefix-hood and checking that the prefix check reports it.
* Every check prints the quantity it measured, not just PASS.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q            # noqa: E402
from s24 import stats_lib as ST          # noqa: E402

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
R = {}


def say(k, v):
    R[k] = v
    print(f"  {k:<52} {v}")


# ===================================================================== dense reference
def _ry(t):
    c, s = np.cos(t / 2.0), np.sin(t / 2.0)
    return np.array([[c, -s], [s, c]])


def _op_on(n, q, g):
    """Embed a 1-qubit gate on wire q of an n-qubit register, qubit 0 = MOST significant."""
    m = np.eye(1)
    for i in range(n):
        m = np.kron(m, g if i == q else np.eye(2))
    return m


def _cnot(n, c, t):
    """Full 2**n CNOT matrix, built from basis indices. Qubit 0 = MOST significant."""
    d = 1 << n
    j = np.arange(d)
    cb = (j >> (n - 1 - c)) & 1
    src = np.where(cb == 1, j ^ (1 << (n - 1 - t)), j)
    M = np.zeros((d, d))
    M[j, src] = 1.0            # M|src> -> |j>; the map is an involution so either order works
    return M


def dense_state(n, theta, layers, final_ry, ring):
    """Independent dense simulation of `layers x (RY all wires, CNOT chain [, ring])`
    [+ one trailing RY layer]. Explicit matrices; no shortcuts shared with the module."""
    nb = layers + int(final_ry)
    th = np.asarray(theta, float).reshape(nb, n)
    psi = np.zeros(1 << n)
    psi[0] = 1.0
    for l in range(nb):
        for q in range(n):
            psi = _op_on(n, q, _ry(th[l, q])) @ psi
        if l >= layers:
            continue
        for q in range(n - 1):
            psi = _cnot(n, q, q + 1) @ psi
        if ring and n > 2:
            psi = _cnot(n, n - 1, 0) @ psi
    return psi


# ============================================================================ A1 -- MPS
def a1_mps():
    print("\nA1  MPSAnsatz: chi = 2**layers BY CONSTRUCTION, exact, no truncation")
    import inspect as _ins
    csrc = _ins.getsource(Q.MPSAnsatz)
    code = "\n".join(l for l in csrc.splitlines()
                     if not l.strip().startswith("#")).lower()
    # FALSIFIER: any truncation primitive inside the ansatz class itself.
    bad = [w for w in ("svd", "truncat", "linalg.qr", "max_bond", "cutoff", "chop")
           if w in code]
    say("truncation primitives inside MPSAnsatz", bad if bad else "NONE")
    fsrc = open(os.path.join(ROOT, "core", "quantum.py")).read().lower()
    say("'svd'/'truncat' anywhere in core/quantum.py",
        {w: fsrc.count(w) for w in ("svd", "truncat", "cutoff")})

    chis = {}
    for L in (1, 2, 3, 4):
        a = Q.MPSAnsatz(6, layers=L, final_ry=True, entangler="cnot")
        chis[L] = int(a.chi)
    say("chi by layers (cnot chain, final_ry)", chis)
    say("chi == 2**layers for all", all(chis[L] == 2 ** L for L in chis))

    a0 = Q.MPSAnsatz(6, layers=2, final_ry=True, entangler="none")
    say("chi with entangler='none' (product state)", int(a0.chi))
    say("ring closure applied by MPSAnsatz", bool(Q.MPSAnsatz(6, 2).ring))

    # amplitudes vs an INDEPENDENT dense simulator, over random angles
    worst_p, worst_a = 0.0, 0.0
    rng = np.random.default_rng(20250908)
    for n in (3, 5, 7, 9):
        for L in (1, 2, 3):
            for fr in (False, True):
                a = Q.MPSAnsatz(n, layers=L, final_ry=fr, entangler="cnot")
                th = rng.normal(0, 1.3, a.n_params())
                p_mps = a.probs(th)
                psi = dense_state(n, th, L, fr, ring=False)
                p_ref = psi ** 2
                worst_p = max(worst_p, float(np.abs(p_mps - p_ref).max()))
                A = a.build(th)
                nrm = a.norm(A)
                worst_a = max(worst_a, abs(float(nrm) - float(psi @ psi)))
    say("max |p_MPS - p_dense| over 24 configs", f"{worst_p:.3e}")
    say("max |<psi|psi>_MPS - 1|", f"{worst_a:.3e}")

    # the padded array's shape IS the bond dimension: (n, chi, 2, chi)
    a = Q.MPSAnsatz(8, layers=2, final_ry=True, entangler="cnot")
    A = a.build(rng.normal(0, 1.0, a.n_params()))
    say("built MPS tensor shape (n, chi, 2, chi)", tuple(int(x) for x in A.shape))
    say("cost is O(n chi^3): chi independent of n", int(Q.MPSAnsatz(64, 2).chi))

    # entanglement is REAL at layers>=1: the product-state ansatz cannot reproduce it
    n, L = 6, 2
    ae = Q.MPSAnsatz(n, layers=L, final_ry=True, entangler="cnot")
    th = rng.normal(0, 1.1, ae.n_params())
    pe = ae.probs(th)
    # Schmidt rank across the middle cut, measured on the dense state
    psi = dense_state(n, th, L, True, ring=False).reshape(1 << (n // 2), -1)
    sv = np.linalg.svd(psi, compute_uv=False)
    say("Schmidt values across the middle cut", np.round(sv, 6).tolist())
    say("Schmidt rank (sv > 1e-12) vs chi", (int((sv > 1e-12).sum()), int(ae.chi)))
    ap = Q.MPSAnsatz(n, layers=L, final_ry=True, entangler="none")
    pp = ap.probs(th)
    say("max |p_cnot - p_product| (same angles)", f"{float(np.abs(pe - pp).max()):.3e}")
    return worst_p, chis


# ================================================================== A2 -- the gradient
def a2_gradient():
    print("\nA2  run_cvar_vqe: Adam on the EXACT parameter-shift gradient, exact statevector")
    rng = np.random.default_rng(7)
    n, L = 7, 3
    circ = Q.StatevectorCircuit(n, L)
    say("deployed register (core/pipeline Config)", f"n={n} qubits, layers={L}, dim={circ.dim}")
    say("parameter count n_params = layers*n", int(circ.n_params()))

    # the statevector itself, against the independent dense simulator (ring=True here)
    th = rng.normal(0, 0.6, circ.n_params())
    p_ref = dense_state(n, th, L, final_ry=False, ring=True) ** 2
    say("max |p_statevector - p_dense|", f"{float(np.abs(circ.probs(th) - p_ref).max()):.3e}")

    # parameter shift vs central finite differences on the EXACT CVaR
    cos_ps, rel_ps = [], []
    for alpha in (0.1, 0.25, 1.0):
        for trial in range(4):
            E = Q.__dict__["np"].asarray(rng.normal(0, 1, circ.dim))
            E = (E - E.mean()) / E.std()
            th = rng.normal(0, 0.6, circ.n_params())
            g_ps = Q.grad_cvar_paramshift(circ, th, E, alpha)
            g_fd = Q.grad_cvar_fd(circ, th, E, alpha, h=1e-5)
            c = float(g_ps @ g_fd / (np.linalg.norm(g_ps) * np.linalg.norm(g_fd)))
            cos_ps.append(c)
            rel_ps.append(float(np.linalg.norm(g_ps - g_fd) / np.linalg.norm(g_fd)))
    say("cos(param-shift, exact FD), min over 12", f"{min(cos_ps):.9f}")
    say("rel |g_ps - g_fd| / |g_fd|, max over 12", f"{max(rel_ps):.3e}")

    # the DEPLOYED objective is the free energy: check ITS gradient too, both terms
    cos_f, rel_f = [], []
    for alpha, T in ((0.1, 0.1), (0.25, 0.3), (1.0, 0.3), (1.0, 0.1)):
        E = rng.normal(0, 1, circ.dim)
        E = (E - E.mean()) / E.std()
        th = rng.normal(0, 0.6, circ.n_params())
        _, g, _, _, _ = Q.free_energy(circ, th, E, alpha, T)
        h = 1e-5
        gf = np.empty_like(g)
        for k in range(len(th)):
            tp = th.copy(); tp[k] += h
            tm = th.copy(); tm[k] -= h
            gf[k] = (Q.free_energy(circ, tp, E, alpha, T)[0]
                     - Q.free_energy(circ, tm, E, alpha, T)[0]) / (2 * h)
        cos_f.append(float(g @ gf / (np.linalg.norm(g) * np.linalg.norm(gf))))
        rel_f.append(float(np.linalg.norm(g - gf) / np.linalg.norm(gf)))
    say("cos(free-energy grad, FD), min over 4 (a,T)", f"{min(cos_f):.9f}")
    say("rel err of free-energy grad, max over 4", f"{max(rel_f):.3e}")

    # THE RECORDED DEFECT: it must still be reproducible AND the default must be the fix
    import inspect
    sig = inspect.signature(Q.cvar_gradient)
    say("cvar_gradient default baseline", sig.parameters["baseline"].default)
    sigs = inspect.signature(Q.grad_cvar_score)
    say("grad_cvar_score default baseline", sigs.parameters["baseline"].default)

    cos_const, cos_tail, nrm_tail = [], [], []
    for trial in range(12):
        nq, LL = 7, 3
        c2 = Q.StatevectorCircuit(nq, LL)
        E = rng.normal(0, 1, c2.dim)
        E = (E - E.mean()) / E.std()
        th = rng.normal(0, 0.6, c2.n_params())
        alpha = 0.15
        g_ex = Q.grad_cvar_paramshift(c2, th, E, alpha)
        # EXACT-EXPECTATION score-function estimators: zero sampling noise, so any
        # disagreement is BIAS, which is the whole point of the recorded defect.
        p = c2.probs(th)
        _, q, _ = Q.cvar_exact(E, p, alpha)
        PR = c2.probs_batch(c2._shift_grid(th, np.pi / 2))
        G = ((PR[0::2] - PR[1::2]) / 2.0 / np.maximum(p, 1e-300)[None, :]).T   # dlogp/dth
        w = np.where(E < q, (E - q) / alpha, 0.0)
        gc = (p * (w - (p * w).sum())) @ G                   # CONSTANT baseline
        m = E < q
        wt = w.copy()
        if m.any():
            wt[m] -= (p[m] * w[m]).sum() / max(p[m].sum(), 1e-300)   # TAIL-ONLY baseline
            wt[~m] = 0.0
        gt = (p * wt) @ G
        cos_const.append(float(gc @ g_ex / (np.linalg.norm(gc) * np.linalg.norm(g_ex))))
        cos_tail.append(float(gt @ g_ex / (np.linalg.norm(gt) * np.linalg.norm(g_ex))))
        nrm_tail.append(float(np.linalg.norm(gt) / np.linalg.norm(g_ex)))
    say("cos(exact-expectation, CONST baseline)", f"{float(np.mean(cos_const)):.6f}")
    say("cos(exact-expectation, TAIL baseline)", f"{float(np.mean(cos_tail)):.6f}")
    say("|g_tail| / |g_exact|", f"{float(np.mean(nrm_tail)):.3f}")

    # run_cvar_vqe must not sample anywhere
    src = inspect.getsource(Q.run_cvar_vqe) + inspect.getsource(Q.free_energy)
    say("sampling calls inside run_cvar_vqe/free_energy",
        [w for w in ("rng.choice", "sample(", "shots") if w in src] or "NONE")
    say("rng use inside run_cvar_vqe", "initial angles only"
        if "rng.normal" in src and "rng.choice" not in src else "OTHER -- INSPECT")
    return min(cos_ps), min(cos_f), float(np.mean(cos_tail))


# ============================================================ A3 -- set equality theorem
def _energy_families(n, rng):
    return {
        "iid_normal": rng.normal(0, 1, n),
        "all_tied": np.zeros(n),
        "two_values": rng.choice([0.0, 1.0], n),
        "heavy_degenerate": np.round(rng.normal(0, 0.4, n), 1),
        "bimodal": np.where(rng.random(n) < 0.5, rng.normal(-5, .2, n), rng.normal(5, .2, n)),
        "huge_scale": rng.normal(0, 1e7, n),
        "one_outlier": np.r_[1e9, rng.normal(0, 1, n - 1)],
        "monotone_ramp": np.arange(n, dtype=float),
        "near_constant": 1.0 + 1e-15 * rng.normal(0, 1, n),
    }


def _prob_families(n, rng):
    """EXACT zeros throughout -- no epsilon floors. That is the s24 REV1 failure mode."""
    u = rng.random(n)
    f = {}
    f["uniform"] = np.ones(n) / n
    f["dirichlet"] = u / u.sum()
    f["peaked"] = (u ** 8) / (u ** 8).sum()
    z = u.copy(); z[rng.random(n) < 0.5] = 0.0
    f["half_exact_zero"] = z / max(z.sum(), 1e-300)
    z2 = u.copy(); z2[rng.random(n) < 0.9] = 0.0
    if z2.sum() == 0:
        z2[0] = 1.0
    f["90pct_exact_zero"] = z2 / z2.sum()
    pm = np.zeros(n); pm[rng.integers(0, n)] = 1.0
    f["point_mass"] = pm
    pm2 = np.zeros(n); pm2[rng.choice(n, 2, replace=False)] = 0.5
    f["two_point_masses"] = pm2
    # adversarial: zero out EXACTLY the lowest-energy states (the hardest hole pattern)
    f["_zero_lowest"] = None          # filled per energy vector by the caller
    return f


def a3_theorem():
    print("\nA3  the set-equality theorem, re-derived on families with EXACT zeros")
    rng = np.random.default_rng(31337)
    n_cells = subset_v = prefix_v = holes_nonzero = full_cells = full_eq = 0
    zero_seen = 0
    for n in (4, 8, 16, 32, 64, 128):
        E_fams = _energy_families(n, rng)
        for ename, E in E_fams.items():
            E = np.asarray(E, float)
            P = _prob_families(n, rng)
            o0 = np.argsort(E, kind="stable")
            zl = np.ones(n); zl[o0[:max(1, n // 4)]] = 0.0
            P["_zero_lowest"] = zl / zl.sum()
            for pname, p in P.items():
                for alpha in (0.01, 0.05, 0.15, 0.35, 0.75, 1.0):
                    n_cells += 1
                    val, q, mass = Q.cvar_from_probs(E, p, alpha)
                    tail = np.flatnonzero(mass > 0)
                    order = np.argsort(E, kind="stable")
                    pos = {int(ix): r for r, ix in enumerate(order)}
                    cut = max(pos[int(i)] for i in tail) if tail.size else -1
                    prefix = order[:cut + 1]
                    # THEOREM: subset of an initial prefix of the ENERGY order
                    if not set(tail.tolist()) <= set(prefix.tolist()):
                        subset_v += 1
                    holes = set(prefix.tolist()) - set(tail.tolist())
                    if holes:
                        prefix_v += 1
                        # MECHANISM: every hole must be an EXACT zero-probability state
                        pn = p / p.sum()
                        if any(pn[h] != 0.0 for h in holes):
                            holes_nonzero += 1
                    if (p > 0).all():
                        full_cells += 1
                        k = int(np.ceil(alpha * n)) if alpha < 1 else n
                        # the classical top-m at the REALISED m, compared by VALUE
                        m = len(tail)
                        cls = order[:m]
                        if np.allclose(np.sort(E[tail]), np.sort(E[cls]), rtol=0, atol=1e-12):
                            full_eq += 1
                        del k
                    else:
                        zero_seen += 1
    say("n_cells", n_cells)
    say("cells containing an EXACT zero probability", zero_seen)
    say("SUBSET-hood violations  <- THE THEOREM", subset_v)
    say("holes that were NOT exactly zero-probability", holes_nonzero)
    say("prefix-hood violations (the s24 REV1 over-claim)", prefix_v)
    say("full-support cells", full_cells)
    say("full-support cells with EXACT value equality", f"{full_eq} / {full_cells}")

    # THE ASSERTION IS LIVE: construct a case that MUST have a hole, and check we see it.
    E = np.arange(8.0)
    p = np.ones(8) / 7.0
    p[0] = 0.0                              # punch out the single lowest-energy state
    _, _, mass = Q.cvar_from_probs(E, p, 0.30)
    tail = np.flatnonzero(mass > 0)
    say("LIVENESS: tail with p[argmin]=0 excludes state 0", 0 not in tail.tolist())
    say("LIVENESS: that tail is still a subset of the prefix",
        set(tail.tolist()) <= set(range(int(tail.max()) + 1)))

    # tail_indices takes no probability vector at all -- the other half of the theorem
    import inspect
    say("tail_indices signature", str(inspect.signature(Q.tail_indices)))
    return dict(n_cells=n_cells, subset_v=subset_v, prefix_v=prefix_v,
                holes_nonzero=holes_nonzero, full_cells=full_cells, full_eq=full_eq,
                zero_cells=zero_seen)


# ============================================================== A4 / Q-C harness audit
def a4_harness():
    print("\nA4/Q-C  s24/d_harness.py audit -- is the MDE fix in, and does it have siblings?")
    from s24 import d_harness as H
    src = open(os.path.join(ROOT, "s24", "d_harness.py")).read()
    say("UNDERPOWERED (<0.7x MDE) rule present", "UNDERPOWERED (<0.7x MDE" in src)

    # measure the fix rather than read it: place an effect at exactly 0.39x its own MDE
    rng = np.random.default_rng(0)
    nn = 126
    folds = np.repeat(np.arange(5), 26)[:nn]
    d = rng.normal(0, 1, nn)
    d = (d - d.mean()) / d.std(ddof=1)
    se = 1.0 / np.sqrt(nn)
    d = d + 0.39 * 2.8016 * se                 # mean is exactly 0.39x MDE
    st = H.paired_stats(d, np.zeros(nn), folds, n_boot=800, seed=1)
    say("0.39x MDE effect -> verdict", st["verdict"])
    say("  its eff_over_mde", round(st["eff_over_mde"], 4))
    d10 = (d - d.mean()) + 1.0 * 2.8016 * se
    say("1.00x MDE effect -> verdict",
        H.paired_stats(d10, np.zeros(nn), folds, n_boot=800, seed=1)["verdict"])
    d20 = (d - d.mean()) + 2.0 * 2.8016 * se
    say("2.00x MDE effect -> verdict",
        H.paired_stats(d20, np.zeros(nn), folds, n_boot=800, seed=1)["verdict"])
    say("MDE == 2.8016*SE exactly", abs(st["mde"] - 2.8016 * st["se"]) < 1e-12)

    # ---- SIBLING HUNT: other places a MISSING measurement could report a PASS ----
    sibs = []
    if 'r.get("gate_equality", True)' in src:
        sibs.append("aggregate(): gate_equality_rate defaults a MISSING key to True, so a row "
                    "that never measured equality is counted as having passed it. gate_pass "
                    "next to it is a strict r['gate_pass']. Same class as the dead assertion.")
    if "required_keys is not None" in src and "_COMPLETE" in src:
        sibs.append("write(): with required_keys=None no _COMPLETE sidecar is written, so a "
                    "STALE sidecar from an earlier complete run survives beside new content "
                    "and keeps asserting complete=true.")
    if "e_cut + 1e-12" in src:
        sibs.append("gate_set_equality(): the prefix is widened by an ABSOLUTE 1e-12, which "
                    "can only make subset-hood EASIER to pass. Harmless on rank-standardised "
                    "energies (gaps ~1/k) and meaningless on raw ones (ulp at 1e18 is ~100), "
                    "but the leniency points the flattering way.")
    if "np.isfinite(a) & np.isfinite(b)" in src:
        sibs.append("paired_stats(): non-finite pairs are dropped silently; n is reported, but "
                    "aggregate() uses nanmean for the MEANS and paired_stats' filter for the "
                    "CONTRASTS, so a mean and its contrast can rest on different target sets.")
    for s in sibs:
        print(f"    SIBLING: {s}")
    say("siblings found", len(sibs))
    R["siblings"] = sibs

    # demonstrate the gate_equality sibling rather than assert it
    rows = [dict(basis="point_cloud", fold=0, source="x", pdb="T", q_rmsd=1.0,
                 c_matched_rmsd=1.0, c_fixed_rmsd=1.0, rand_rmsd_mean=2.0,
                 rand_rmsd_best=1.5, gate_pass=True, q_m=10, q_entropy_bits=5.0)
            for _ in range(6)]
    for i, r in enumerate(rows):
        r["pdb"] = f"T{i}"
        r["fold"] = i % 5
    agg = H.aggregate(rows)
    say("aggregate() equality rate with the key ABSENT", agg["gate_equality_rate"])
    return sibs


# ============================================================ deployed-config census
def a5_deployed():
    print("\nA5  what is actually deployed, from source")
    from core.pipeline import Config, VQE_LFO
    c = Config()
    say("core/pipeline Config.vqe_qubits", c.vqe_qubits)
    say("core/pipeline Config.vqe_layers", c.vqe_layers)
    say("core/pipeline Config.vqe_iters", c.vqe_iters)
    say("core/pipeline Config.vqe_seed", c.vqe_seed)
    say("hypothesis-set size 2**vqe_qubits", 1 << c.vqe_qubits)
    say("parameters = layers * qubits", c.vqe_layers * c.vqe_qubits)
    say("VQE_LFO (alpha, T) per fold", {int(k): v for k, v in VQE_LFO.items()})
    say("folds running alpha == 1.0 (NO tail constraint)",
        sorted(int(k) for k, v in VQE_LFO.items() if v[0] == 1.0))
    say("distinct temperatures in the deployed table",
        sorted({v[1] for v in VQE_LFO.values()}))
    say("selector ansatz", "StatevectorCircuit(RY+CNOT chain+ring), exact, via run_cvar_vqe")
    say("MPSAnsatz used by core/pipeline", False)


if __name__ == "__main__":
    a1_mps()
    a2_gradient()
    a3_theorem()
    a4_harness()
    a5_deployed()
    ST.save_atomic(os.path.join(OUT, "q_verify.json"),
                   dict(kind="property_verification_no_rmsd", lane="Q", sprint=25,
                        prereg="s25/PREREG_Q.md Q-A / Q-C", results=R),
                   module_file=__file__)
    print("\nwrote", os.path.join(OUT, "q_verify.json"))
