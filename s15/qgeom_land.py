"""SPRINT 15 / QGEOM -- PART E: THE LANDSCAPE.  Exact Hessians at three points of a run.

For the expectation-value cost ``C(theta) = sum_x E(x) p(x)``, ``p = psi^2`` with real
amplitudes, the exact Hessian follows from the same shift identity used everywhere in this
workstream:

    d/dt RY(t)   = (1/2) RY(t + pi)
    d2/dt2 RY(t) = (1/4) RY(t + 2pi) = -(1/4) RY(t)          [RY has period 4pi]

so, writing `psi_{i+j}` for the state with BOTH theta_i and theta_j advanced by pi,

    d_i d_j psi = (1/4) psi_{i+j}   (i != j),     d_i^2 psi = -(1/4) psi
    H_ij = 2 sum_x E(x) [ (d_i d_j psi)(x) psi(x) + (d_i psi)(x) (d_j psi)(x) ]

Exact to machine precision, `P(P-1)/2 + P` statevector builds, no finite differences.
Cross-checked against central differences of the exact gradient.

WHAT IS MEASURED at initialisation, mid-optimisation and at the final point:
leading and smallest eigenvalues, the number of NEGATIVE directions, effective rank,
condition number, the flat-direction count, the gradient norm, and the alignment of the
gradient with the Hessian's eigenbasis.

THE DISCIPLINE THE BRIEF INSISTS ON: *do not confuse optimisation convergence with reaching
a meaningful basin.*  So every row prints the curvature evidence that a point IS a minimum
next to the structural quality of what that minimum contains.  A clean positive-definite
minimum holding a 5 A structure is the finding, not a failure of the optimiser.

    python -m s15.qgeom_land
"""
from __future__ import annotations

import time

import numpy as np

from core.quantum import Adam, cvar_exact
from s15 import qgeom_lib as G
from s15.qgeom_qng import sub_objective

TAG = "land"


def hessian_expectation(circ, theta, E):
    """Exact Hessian of the expectation value.  O(P^2) statevector builds."""
    th = np.asarray(theta, float)
    P = th.size
    psi = circ.state(th)
    D = G.deriv_states(circ, th)
    H = 2.0 * (D * E[None, :]) @ D.T                       # the (d_i psi)(d_j psi) term
    Epsi = E * psi
    for i in range(P):
        H[i, i] += 2.0 * (-0.25) * float(psi @ Epsi)
        for j in range(i + 1, P):
            t = th.copy()
            t[i] += np.pi
            t[j] += np.pi
            v = 0.25 * circ.state(t)
            c = 2.0 * float(v @ Epsi)
            H[i, j] += c
            H[j, i] += c
    return (H + H.T) / 2.0


def hess_stats(H, grad=None, tol=1e-9):
    ev = np.linalg.eigvalsh(H)
    P = H.shape[0]
    scale = max(np.abs(ev).max(), 1e-30)
    flat = int(np.sum(np.abs(ev) <= tol * scale))
    neg = int(np.sum(ev < -tol * scale))
    pos = np.abs(ev)[np.abs(ev) > tol * scale]
    out = {"n_params": P, "eig_max": float(ev.max()), "eig_min": float(ev.min()),
           "n_negative": neg, "n_flat": flat,
           "frac_negative": float(neg / P), "frac_flat": float(flat / P),
           "abs_cond": float(pos.max() / pos.min()) if pos.size else float("nan"),
           "trace": float(np.trace(H)),
           "eff_rank": float(np.exp(-(lambda s: (s * np.log(s)).sum())(
               np.abs(ev)[np.abs(ev) > 0] / np.abs(ev).sum()))),
           "eigs": [float(x) for x in ev]}
    if grad is not None:
        w, V = np.linalg.eigh(H)
        comp = (V.T @ grad) ** 2
        tot = max(comp.sum(), 1e-300)
        out["grad_norm"] = float(np.linalg.norm(grad))
        out["grad_share_negative_curv"] = float(comp[w < -tol * scale].sum() / tot)
        out["grad_share_flat"] = float(comp[np.abs(w) <= tol * scale].sum() / tot)
    return out


def verify_hessian(n=8, pat="ring", L=2, seed=0, h=1e-5):
    """Independent check against central differences of the EXACT gradient."""
    c = G.FlexCircuit(n, L, pat, 2)
    E = G.V.uniformise(np.random.default_rng(0).standard_normal(c.dim))
    th = G.random_theta(c, seed)
    H = hessian_expectation(c, th, E)

    def grad(t):
        psi = c.state(t)
        return G.grad_adjoint(c, t, 2.0 * psi * E)

    Hfd = np.empty_like(H)
    for j in range(th.size):
        tp = th.copy(); tp[j] += h
        tm = th.copy(); tm[j] -= h
        Hfd[:, j] = (grad(tp) - grad(tm)) / (2 * h)
    Hfd = (Hfd + Hfd.T) / 2
    return float(np.abs(H - Hfd).max() / max(np.abs(H).max(), 1e-30))


def run(pdbs=("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5"),
        residues=(1, 2, 3, 4, 5, 6), patterns=("ring", "block"), L=2,
        iters=200, seeds=(0, 1), lr=0.10):
    print("=" * 118)
    print("E1. EXACT HESSIAN AT INITIALISATION, MID-OPTIMISATION AND AT THE FINAL POINT")
    print("=" * 118)
    v = verify_hessian()
    print(f"  verification: exact Hessian vs central differences of the exact gradient, "
          f"max rel err {v:.3e}")
    G.ck(TAG, "verify_rel_err", v)
    prev = G.ck_load(TAG)
    out = prev.get("E1", {})
    print()
    print(f"{'cell':26s} {'point':8s} {'|grad|':>9s} {'lam_max':>9s} {'lam_min':>10s} "
          f"{'#neg':>5s} {'#flat':>6s} {'|cond|':>10s} {'effrank':>8s} "
          f"{'g in neg':>9s} {'cost':>8s} {'modeRMSD':>9s}")
    for pdb in pdbs:
        e = G.V.Enum(pdb)
        base = np.random.default_rng(0).integers(0, e.k, e.n)
        E = G.V.uniformise(sub_objective(e, residues, base, "legacy"))
        R = sub_objective(e, residues, base, "rmsd")        # ORACLE, post-hoc
        for pat in patterns:
            c = G.FlexCircuit(2 * len(residues), L, pat, 2)
            for s in seeds:
                key = f"{pdb}|{pat}|s{s}"
                if key in out:
                    continue
                t0 = time.time()
                th = G.random_theta(c, s)
                opt = Adam(c.n_params(), lr=lr)
                snaps = {}
                for it in range(iters + 1):
                    psi = c.state(th)
                    p = psi ** 2
                    p = p / p.sum()
                    gr = G.grad_adjoint(c, th, 2.0 * psi * E)
                    if it in (0, iters // 2, iters):
                        nm = {0: "init", iters // 2: "mid", iters: "final"}[it]
                        H = hessian_expectation(c, th, E)
                        st = hess_stats(H, gr)
                        st["cost"] = float(p @ E)
                        st["mode_rmsd"] = float(R[int(np.argmax(p))])
                        st["metric"] = G.spec_stats(G.fs_metric(c, th))
                        st.pop("eigs", None)
                        st["metric"].pop("eigs", None)
                        snaps[nm] = st
                    if it < iters:
                        th = opt.step(th, gr)
                out[key] = snaps
                for nm in ("init", "mid", "final"):
                    st = snaps[nm]
                    print(f"{key:26s} {nm:8s} {st['grad_norm']:9.5f} "
                          f"{st['eig_max']:9.4f} {st['eig_min']:10.4f} "
                          f"{st['n_negative']:5d} {st['n_flat']:6d} "
                          f"{st['abs_cond']:10.2f} {st['eff_rank']:8.2f} "
                          f"{st['grad_share_negative_curv']:9.4f} {st['cost']:8.4f} "
                          f"{st['mode_rmsd']:9.3f}")
                print(f"{'':26s} ({time.time()-t0:.0f}s)")
                G.ck(TAG, "E1", out)
    summarise(out)
    return out


def summarise(out):
    print()
    print("=" * 118)
    print("E2. SUMMARY -- convergence is NOT the same as a meaningful basin")
    print("=" * 118)
    res = {}
    for nm in ("init", "mid", "final"):
        rows = [v[nm] for v in out.values() if nm in v]
        if not rows:
            continue
        f = lambda k: float(np.mean([r[k] for r in rows]))
        res[nm] = {k: f(k) for k in ("grad_norm", "eig_max", "eig_min", "n_negative",
                                     "n_flat", "abs_cond", "eff_rank",
                                     "grad_share_negative_curv", "cost", "mode_rmsd")}
        res[nm]["frac_cells_PSD"] = float(np.mean([r["n_negative"] == 0 for r in rows]))
        res[nm]["n_cells"] = len(rows)
        print(f"  {nm:6s} n={len(rows):3d}  |grad| {f('grad_norm'):.5f}  "
              f"lam_max {f('eig_max'):+.4f}  lam_min {f('eig_min'):+.4f}  "
              f"#neg {f('n_negative'):.2f}  #flat {f('n_flat'):.2f}  "
              f"PSD in {res[nm]['frac_cells_PSD']:.1%} of cells  "
              f"cost {f('cost'):.4f}  mode RMSD {f('mode_rmsd'):.3f} A")
    G.ck(TAG, "E2_summary", res)
    return res


def main():
    G.wait_mem(0.8, "qgeom_land")
    G.ck_load(TAG)
    run()
    print("\nwritten -> s15/results/qgeom_land.json")


if __name__ == "__main__":
    main()
