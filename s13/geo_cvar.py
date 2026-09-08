"""SPRINT 13 GEO -- EXPERIMENT 6: CVaR as a research variable, not an ingredient.

alpha in {1.0, 0.5, 0.25, 0.1, 0.05, 0.01} x {Legacy, AMBER, AMBER-SOFT}, everything else
matched (representation, ansatz, theta0, optimiser, steps, learning rate, seeds).

WHAT IS MEASURED, AND THE THREE CONTROLS THAT MAKE IT MEAN SOMETHING
====================================================================
Measured per arm: gradient variance at random theta, the optimisation trace, final CVaR,
state entropy, ``p_top``, the probability mass the state puts on the lowest-energy 1% of the
register ("low-energy tail mass"), and -- separately, never merged with the objective --
the ORACLE CA-RMSD of three readouts.

CONTROL 1 -- **matched-entropy Boltzmann.**  s12 measured that CVaR's breadth is
*indifference rather than preference*: its distribution lost to a matched-entropy Boltzmann
distribution on 20 of 20 arms.  So every VQE readout here is scored against
``p ∝ exp(-E/T')`` with ``T'`` bisected so that ``H(Boltzmann) == H(VQE state)``, read out
identically.  Beating a one-line softmax is the bar.

CONTROL 2 -- **estimability, not just trainability.**  Evaluation-budget parity is NOT shot
parity once alpha differs: at alpha = 0.05 only 5 % of the drawn samples enter the CVaR, so
the estimator's effective sample size is ``alpha * shots`` -- 20x smaller at the same
evaluation count.  The empirical sd of the sampled CVaR estimate at fixed shots is measured
per alpha, so "small alpha destroys trainability" and "small alpha destroys estimability"
are separated.  (Qiu et al. 2026: tilting moves the bottleneck from trainability to
estimability and does NOT remove barren plateaus.  **Nothing here describes CVaR as a
barren-plateau mitigation.**)

CONTROL 3 -- **AMBER-SOFT.**  A monotone (rank-preserving) compression of AMBER's upper
tail.  If an alpha effect on AMBER disappears under AMBER-SOFT it was dynamic range; if it
survives it was the ranking.

    python -m s13.geo_cvar
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from s12 import instrument as I                        # noqa: E402
from core.quantum import Adam, cvar_exact              # noqa: E402

ALPHAS = (1.0, 0.5, 0.25, 0.1, 0.05, 0.01)
LAYERS = 3
STEPS = 200
LR = 0.10
SEEDS = (0, 1)
SHOTS = 1024
CELLS = [(p, 5, 4) for p in ("1A13", "1A1P", "2BFI")] + [(p, 6, 4) for p in ("1A13", "1A1P", "2BFI")]


def ca_coords(pdb, L, k):
    """(N, L, 3) CA traces of every register state.  ORACLE-free (structure, not score)."""
    t = G.peptide(pdb, L)
    rep = G.make_rep(t, k)
    S = G.all_states(L, k)
    rows = np.arange(L)
    PHI = np.asarray(rep._phi, float)[rows[None, :], S]
    PSI = np.asarray(rep._psi, float)[rows[None, :], S]
    out = np.empty((len(S), L, 3))
    for a in range(0, len(S), 4096):
        out[a:a + 4096] = I.build_ca(PHI[a:a + 4096], PSI[a:a + 4096])
    return out, t["nat_ca"]


def readouts(p, R, W, nat, top=256):
    """ORACLE EVALUATION.  Three readouts of a distribution over register states."""
    idx = np.argsort(p)[::-1][:top]
    q = p[idx] / p[idx].sum()
    avg = (W[idx] * q[:, None, None]).sum(0)
    return {"rmsd_argmax_p": float(R[int(np.argmax(p))]),
            "rmsd_best_in_top16": float(R[np.argsort(p)[::-1][:16]].min()),
            "rmsd_p_weighted_top256": float(I.ca_rmsd(avg, nat))}


def gibbs_matched(E, H_target, lo=1e-6, hi=1e12, iters=80):
    """Boltzmann p ∝ exp(-E/T) with T bisected so its entropy matches `H_target` (bits)."""
    E = np.asarray(E, float)
    Em = E - E.min()

    def ent(T):
        z = -Em / max(T, 1e-300)
        z -= z.max()
        p = np.exp(z); p /= p.sum()
        return G.entropy_bits(p), p

    hH, hp = ent(hi)
    if H_target >= hH:
        return hp, hi
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        h, p = ent(mid)
        if h < H_target:
            lo = mid
        else:
            hi = mid
    return ent(np.sqrt(lo * hi))[1], float(np.sqrt(lo * hi))


def grad_var_at_random(E, circ, alpha, n_theta=32, seed=0):
    rng = np.random.default_rng(seed)
    Gd = []
    for _ in range(n_theta):
        th = rng.uniform(-np.pi, np.pi, circ.n_params())
        Gd.append(G.cvar_grad(E, circ, th, alpha))
    Gd = np.array(Gd)
    iqr = G.scales(E)["iqr"] or 1.0
    return {"var_mean": float(Gd.var(0).mean()),
            # IQR-normalised: removes the energy model's units so Legacy and AMBER are
            # comparable.  `var_mean_iqr_alpha` additionally removes the trivial 1/alpha in
            # the envelope-theorem weight ``(q-E)_+/alpha``, so an alpha effect that
            # survives it is a landscape effect rather than a rescaling.
            "var_mean_iqr": float(Gd.var(0).mean() / iqr ** 2),
            "var_mean_iqr_alpha": float(Gd.var(0).mean() * alpha ** 2 / iqr ** 2),
            "grad_norm_mean": float(np.linalg.norm(Gd, axis=1).mean()),
            "grad_norm_mean_iqr": float(np.linalg.norm(Gd, axis=1).mean() / iqr)}


def estimability(E, circ, theta, alpha, shots=SHOTS, reps=64, seed=0):
    """Shot noise of the SAMPLED CVaR at fixed shots -- the estimability axis."""
    rng = np.random.default_rng(seed)
    p = circ.probs(theta)
    exact = float(cvar_exact(E, p, alpha)[0])
    vals = []
    for _ in range(reps):
        x = rng.choice(len(p), size=shots, p=p)
        e = np.sort(E[x])
        m = max(1, int(round(alpha * shots)))
        vals.append(float(e[:m].mean()))
    vals = np.array(vals)
    # The scale must not be `|CVaR - min(E)|`: a converged state drives that to zero and the
    # ratio to +inf, which is an artefact of convergence rather than a shot-noise statement.
    # IQR(E) is a fixed property of the register, so `relative_sd` is comparable across
    # alphas, models and cells.
    iqr = G.scales(E)["iqr"]
    return {"shots": shots, "effective_sample_size": float(alpha * shots),
            "cvar_exact": exact, "cvar_sampled_mean": float(vals.mean()),
            "cvar_sampled_sd": float(vals.std()),
            "cvar_sampled_bias": float(vals.mean() - exact),
            "relative_sd": float(vals.std() / max(iqr, 1e-300)),
            "relative_bias": float((vals.mean() - exact) / max(iqr, 1e-300))}


def run_cell(pdb, L, k):
    n = int(round(np.log2(k))) * L
    circ = G.circuit(n, layers=LAYERS)
    P = circ.n_params()
    R = G.rmsd_table(pdb, L, k)
    W, nat = ca_coords(pdb, L, k)
    lowmask = R * 0
    rows = []
    for v in G.VARIANTS:
        E = G.variant_table(pdb, L, k, v)
        thr = float(np.quantile(E, 0.01))
        low = E <= thr
        for a in ALPHAS:
            gv = grad_var_at_random(E, circ, a)
            for s in SEEDS:
                th0 = G.init_theta(P, s)
                th = th0.copy()
                opt = Adam(P, lr=LR)
                trace = []
                for _ in range(STEPS):
                    trace.append(G.cvar_value(E, circ, th, a))
                    th = opt.step(th, G.cvar_grad(E, circ, th, a))
                p = circ.probs(th)
                Hb = G.entropy_bits(p)
                pg, T = gibbs_matched(E, Hb)
                rows.append({
                    "pdb": pdb, "L": L, "k": k, "n_qubits": n, "layers": LAYERS,
                    "model": v, "alpha": a, "seed": s, "steps": STEPS, "lr": LR,
                    "grad_var_random_theta": gv,
                    "cvar_start": trace[0],
                    "cvar_end": float(G.cvar_value(E, circ, th, a)),
                    "cvar_drop": float(trace[0] - G.cvar_value(E, circ, th, a)),
                    "converged_frac_of_final": float(
                        (trace[0] - trace[len(trace) // 2])
                        / max(trace[0] - G.cvar_value(E, circ, th, a), 1e-300)),
                    "entropy_bits": Hb, "entropy_max_bits": float(n),
                    "p_top": float(p.max()),
                    "low_energy_tail_mass": float(p[low].sum()),
                    "low_energy_tail_mass_uniform": float(low.mean()),
                    "estimability": estimability(E, circ, th, a),
                    "ORACLE_vqe": readouts(p, R, W, nat),
                    "ORACLE_gibbs_matched_entropy": readouts(pg, R, W, nat),
                    "gibbs_T": T,
                    "ORACLE_rmsd_energy_argmin": float(R[int(np.argmin(E))]),
                    "ORACLE_rmsd_pool_best": float(R.min()),
                })
        print("  %-5s L=%-2d k=%d %-11s done" % (pdb, L, k, v), flush=True)
    return rows


def main():
    rows = []
    for pdb, L, k in CELLS:
        if not os.path.exists(os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")):
            print(f"  skip {pdb} L={L} k={k}")
            continue
        G.gate(tries=1)
        t0 = time.time()
        rows += run_cell(pdb, L, k)
        print("  cell %s L=%d k=%d in %.0fs" % (pdb, L, k, time.time() - t0), flush=True)
        G.write("geo_cvar", {"what": "CVaR alpha sweep with matched-entropy Boltzmann and "
                                     "estimability controls", "rows": rows})
    G.write("geo_cvar", {"what": "CVaR alpha sweep with matched-entropy Boltzmann and "
                                 "estimability controls", "rows": rows})


if __name__ == "__main__":
    main()
