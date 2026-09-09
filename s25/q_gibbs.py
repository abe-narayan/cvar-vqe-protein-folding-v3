"""S25 / LANE Q -- HOW CLOSE IS p_theta TO ITS OWN ANALYTIC OPTIMUM? A DIRECT MEASUREMENT.

PROPERTY MEASUREMENT. Reads no RMSD, no native. Registered FALSIFIER, not a fork list.

WHY THIS EXISTS
===============
A draft of `s25/QUANTUM.md` s0(4) said the variational state "attains its own analytic
optimum". What Q-B actually measured was that the circuit arm and the exact-Boltzmann arm have
indistinguishable RMSD at n=126 (0.24-0.42x MDE). **That is a statement about a downstream
endpoint, and an underpowered one; the sentence is a statement about DISTRIBUTIONS.** Claiming
a mechanism from an endpoint is the exact error that produced this project's L2. So the
distributional claim is measured here directly, and the wording is earned or cut on the result.

THE IDENTITY THAT MAKES IT EXACT
================================
For `F(p) = E_p[E] - T H(p)` with H the Shannon entropy in NATS, the unconstrained minimiser
over the simplex is the Gibbs state `p*(x) = exp(-E(x)/T)/Z`, and for ANY p

    F(p) - F(p*)  =  T * KL(p || p*)      exactly.

So the free-energy suboptimality IS the KL divergence up to the factor T. No estimation and no
sampling: both come out of the same exact statevector the deployed driver uses. The identity is
ASSERTED at runtime on every row, so it is a cross-check on both computations rather than an
assumption.

`cvar_exact` at alpha = 1 returns the full mean of E under p, so the deployed objective at
alpha = 1 IS the F above. The three folds of `VQE_LFO` that run alpha = 1.0 are therefore
solving a Gibbs variational problem with a known closed-form answer, and this module measures
how well 21 parameters solve it.

THE SCALE, WHICH IS THE POINT OF THE COMPARATORS
================================================
A divergence is meaningless without units, so every trained state is reported beside
deliberately imperfect states on the same axis: an untrained random theta, a 5-step and a
15-step Adam state (the same run, stopped early), the uniform distribution, and the collapsed
point mass at the argmin. The reader can then see where 50 Adam steps sits on that ladder.

THE MANDATORY CONTROL
=====================
`concentration-is-wrong-when-discrimination-binds` requires that a variational arm be
differenced against **best-of-N from the SAME untrained distribution**, never against an
initialisation MEAN. That control is run here on the objective (N=200, identical theta law).

REGISTERED FALSIFIER: KL(p_theta || p*) > 0.1 nat makes "attains its own analytic optimum"
unearned and the wording must be cut from the specification.

    >>> RESULT: THE FALSIFIER FIRED. KL = 0.902 nats at the deployed temperature. The wording
    >>> was cut. What replaced it is in this module's output and is a different, weaker, and
    >>> more interesting statement.
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

from core import quantum as Q                  # noqa: E402
from core.pipeline import Config, VQE_LFO, _zrank   # noqa: E402
from s24 import stats_lib as ST                # noqa: E402

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

cfg = Config()
NQ, LAY, ITERS, SEED = cfg.vqe_qubits, cfg.vqe_layers, cfg.vqe_iters, cfg.vqe_seed
DIM = 1 << NQ
_r = np.arange(1, DIM + 1, dtype=float)
E = (_r - _r.mean()) / _r.std()          # the deployed energy vector's shape (see s3.1)
R = {"n_qubits": NQ, "layers": LAY, "iters": ITERS, "seed": SEED, "dim": DIM}


def gibbs(T):
    w = -E / T
    w = w - w.max()
    g = np.exp(w)
    return g / g.sum()


def F(p, T):
    p = np.asarray(p, float)
    p = p / p.sum()
    return float((p * E).sum() + T * (p * np.log(np.maximum(p, 1e-300))).sum())


def stats(p, T):
    p = np.asarray(p, float)
    p = p / p.sum()
    g = gibbs(T)
    kl = float((p * np.log(np.maximum(p, 1e-300) / np.maximum(g, 1e-300))).sum())
    return dict(F=F(p, T), kl_nats=kl, kl_bits=kl / math.log(2),
                tv=float(0.5 * np.abs(p - g).sum()),
                entropy_bits=float(-(p * np.log2(np.maximum(p, 1e-300))).sum()))


def trained(T, iters):
    """The deployed driver, stopped at `iters` Adam steps. Same code path, same seed."""
    p, _, _, _ = Q.run_cvar_vqe(E, 1.0, T, n=NQ, layers=LAY, iters=int(iters), seed=SEED)
    return np.asarray(p, float)


if __name__ == "__main__":
    circ = Q.StatevectorCircuit(NQ, LAY)
    print(f"deployed call: run_cvar_vqe(E, alpha, T, n={NQ}, layers={LAY}, iters={ITERS}, "
          f"seed={SEED});  P = {circ.n_params()} parameters, dim = {DIM}")
    print("alpha = 1 ONLY -- that is where a closed-form optimum exists.\n")

    print("1. THE TRAINED STATE AGAINST ITS OWN ANALYTIC OPTIMUM, WITH A SCALE")
    rows = []
    for T in (0.1, 0.3, 1.0):
        g = gibbs(T)
        rng = np.random.default_rng(SEED)
        ladder = {
            "random theta (untrained)": circ.probs(rng.normal(0.0, 0.6, circ.n_params())),
            "5 Adam steps": trained(T, 5),
            "15 Adam steps": trained(T, 15),
            f"{ITERS} Adam steps (DEPLOYED)": trained(T, ITERS),
            "uniform over 2**n": np.ones(DIM) / DIM,
            "point mass at argmin (the collapse)": np.eye(DIM)[int(np.argmin(E))],
            "the Gibbs optimum itself": g,
        }
        dep = any(v == (1.0, T) for v in VQE_LFO.values())
        print(f"\n  T = {T}{'   <- DEPLOYED on folds ' + str(sorted(k for k, v in VQE_LFO.items() if v == (1.0, T))) if dep else ''}")
        print(f"    {'state':<38} {'F':>10} {'KL nats':>10} {'KL bits':>9} {'TV':>8} {'H bits':>8}")
        cell = {}
        for name, p in ladder.items():
            s = stats(p, T)
            cell[name] = s
            kl = "inf" if not np.isfinite(s["kl_nats"]) else f"{s['kl_nats']:10.6f}"
            klb = "inf" if not np.isfinite(s["kl_bits"]) else f"{s['kl_bits']:9.5f}"
            print(f"    {name:<38} {s['F']:10.6f} {kl:>10} {klb:>9} {s['tv']:8.5f} "
                  f"{s['entropy_bits']:8.4f}")
        # the identity is a CHECK, not an assumption
        d = cell[f"{ITERS} Adam steps (DEPLOYED)"]
        resid = abs((d["F"] - cell["the Gibbs optimum itself"]["F"]) - T * d["kl_nats"])
        assert resid < 1e-9, resid
        rows.append(dict(T=T, deployed_at_this_T=dep, identity_residual=float(resid),
                         ladder={k: v for k, v in cell.items()}))
    R["divergence_ladder"] = rows
    print("\n  identity  F(p) - F(p*) == T * KL(p || p*)  asserted < 1e-9 on every trained row")

    print("\n2. THE MANDATORY CONTROL: best-of-N from the UNTRAINED circuit, N = 200")
    print("   (`concentration-is-wrong-when-discrimination-binds`: never difference a")
    print("    variational arm against an initialisation MEAN.)")
    print(f"    {'T':>5} {'F_trained':>11} {'F_init_mean':>12} {'F_init_best':>12} "
          f"{'F_gibbs':>10} {'gap closed':>11} {'beats best-of-N':>16}")
    ctrl = []
    for T in (0.1, 0.3, 1.0):
        Ft = F(trained(T, ITERS), T)
        Fg = F(gibbs(T), T)
        rng = np.random.default_rng(0)
        Fs = np.array([F(circ.probs(rng.normal(0.0, 0.6, circ.n_params())), T)
                       for _ in range(200)])
        closed = (Fs.mean() - Ft) / (Fs.mean() - Fg)
        beats = bool(Ft < Fs.min())
        print(f"    {T:5.2f} {Ft:11.6f} {Fs.mean():12.6f} {Fs.min():12.6f} {Fg:10.6f} "
              f"{closed:11.3f} {str(beats):>16}")
        ctrl.append(dict(T=T, F_trained=Ft, F_init_mean=float(Fs.mean()),
                         F_init_best=float(Fs.min()), F_gibbs=Fg,
                         frac_gap_closed=float(closed), beats_best_of_N=beats,
                         margin=float(Fs.min() - Ft), N=200))
    R["training_control"] = ctrl

    print("\n3. IS THE DEPLOYED HAMILTONIAN'S SPECTRUM TARGET-INDEPENDENT?")
    print("   `top = argsort(sc)` (core/pipeline.py:768), so `sc[top[:128]]` is ascending and")
    print("   `_zrank` of an ascending vector is the standardised ranks 1..128 -- EXCEPT where")
    print("   the score TIES, because rankdata averages ties. Measured on 8 real targets:")
    from s12 import instrument as I
    ref = E
    devs = []
    for pdb in [t["pdb"] for t in I.targets()][:8]:
        t = {x["pdb"]: x for x in I.targets()}[pdb]
        u = I.load_univ(pdb)
        pi = I.pool_idx(u, k=500)
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(int(t["n"]))
        D = I.pair_dists(u["W"][pi], i, j)
        sc = np.asarray(I.shipped_score(dg, D.astype(np.float32).astype(float)), float)
        o = np.argsort(sc, kind="stable")[:DIM]
        d = float(np.abs(_zrank(sc[o]) - ref).max())
        devs.append(dict(pdb=pdb, max_abs_dev=d))
        print(f"     {pdb}   max |E - standardised ranks| = {d:.3e}")
    mx = max(x["max_abs_dev"] for x in devs)
    rng_E = float(E.max() - E.min())
    print(f"   worst deviation {mx:.3e} against an E range of {rng_E:.4f} "
          f"= {mx / rng_E:.2%} of the range; exactly identical on "
          f"{sum(1 for x in devs if x['max_abs_dev'] == 0.0)}/8 targets")
    print("   => the spectrum of H is target-independent TO WITHIN TIE-AVERAGING, not exactly.")
    R["spectrum_target_independence"] = dict(targets=devs, worst=mx, E_range=rng_E,
                                             worst_frac_of_range=mx / rng_E,
                                             n_exactly_identical=int(sum(
                                                 1 for x in devs if x["max_abs_dev"] == 0.0)),
                                             n_checked=len(devs))

    print("\n4. FALSIFIER VERDICT")
    for r_ in rows:
        d = r_["ladder"][f"{ITERS} Adam steps (DEPLOYED)"]
        v = "EARNED" if d["kl_nats"] <= 0.1 else "UNEARNED -- wording CUT from QUANTUM.md"
        print(f"    T={r_['T']:.2f}  KL = {d['kl_nats']:.6f} nats, TV = {d['tv']:.4f}  -> {v}")
    R["falsifier"] = dict(threshold_nats=0.1, fired_against_the_draft=True)

    ST.save_atomic(os.path.join(OUT, "q_gibbs.json"),
                   dict(kind="property measurement, distributional, no RMSD", lane="Q",
                        sprint=25, prereg="s25/PREREG_Q.md Q-A (property, falsifier registered)",
                        results=R), module_file=__file__)
    print("\nwrote", os.path.join(OUT, "q_gibbs.json"))
