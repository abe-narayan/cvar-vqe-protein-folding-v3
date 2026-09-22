"""L1 check: is the DEPLOYED CVaR free energy F(p) = CVaR_alpha(E;p) - T*H(p) exactly
solvable classically in closed form?

Derivation (Rockafellar-Uryasev 1999/2000 variational form + Gibbs variational principle):

  CVaR_alpha^low(E;p) = max_s { s - (1/alpha) sum_i p_i (s - E_i)_+ }        (R-U, lower tail)

  F(p) = CVaR^low - T H(p) = max_s { s + sum_i p_i [ T log p_i - (s-E_i)_+/alpha ] }

  affine-in-p inside a max over s, plus T*sum p log p  =>  F is CONVEX in p
  (strictly convex for T>0).  Concave in s.  Sion's minimax => swap:

  min_p F = max_s { s - T log sum_i exp( (s - E_i)_+ / (alpha T) ) }         (1-D concave)

  and the minimiser is the HINGED GIBBS distribution

  p*_i  propto  exp( (s* - E_i)_+ / (alpha T) )

i.e. UNIFORM on every candidate at or above the VaR level s*, exponentially tilted below it.
The whole 2**n-dimensional optimisation is determined by ONE scalar s*.

This script checks that claim against core.quantum's own cvar_exact and against the shipped
run_cvar_vqe.
"""
import sys, os, time
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))))

from core.quantum import cvar_exact, StatevectorCircuit, free_energy, run_cvar_vqe


def F_of_p(E, p, alpha, T):
    v, q, _ = cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-300))
    H = float(-(p * lp).sum())
    return float(v - T * H), v, H


def dual_g(E, s, alpha, T):
    """g(s) = s - T logsumexp( (s-E)_+/(alpha T) ).  Concave in s."""
    u = np.maximum(s - E, 0.0) / (alpha * T)
    m = u.max()
    return s - T * (m + np.log(np.exp(u - m).sum()))


def closed_form(E, alpha, T, iters=200):
    """Golden-section on the concave 1-D dual; returns (p_star, F_star, s_star)."""
    lo, hi = float(E.min()) - 1e-9, float(E.max()) + 1e-9
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = dual_g(E, c, alpha, T), dual_g(E, d, alpha, T)
    for _ in range(iters):
        if fc < fd:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = dual_g(E, d, alpha, T)
        else:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = dual_g(E, c, alpha, T)
    s = 0.5 * (a + b)
    u = np.maximum(s - E, 0.0) / (alpha * T)
    u -= u.max()
    p = np.exp(u)
    p /= p.sum()
    return p, dual_g(E, s, alpha, T), s


def projected_gd(E, alpha, T, iters=4000, lr=0.02, seed=0):
    """Independent classical reference: mirror descent on the simplex."""
    rng = np.random.default_rng(seed)
    D = len(E)
    w = rng.normal(0, 0.1, D)
    for _ in range(iters):
        w -= w.max()
        p = np.exp(w); p /= p.sum()
        _, _, dCV = cvar_exact(E, p, alpha)
        g = dCV + T * (np.log(np.maximum(p, 1e-300)) + 1.0)
        # mirror step in log space
        gp = p * (g - float(p @ g))
        w -= lr * gp / (p + 1e-12)
    w -= w.max()
    p = np.exp(w); p /= p.sum()
    return p


def main():
    rows = []
    for n, layers in ((7, 3), (9, 3)):
        D = 1 << n
        for alpha in (0.1, 0.25, 1.0):
            for T in (0.1, 0.05):
                rng = np.random.default_rng(11)
                # the deployed E is a rank-standardised score vector: z of ranks
                r = np.arange(D, dtype=float)
                E = (r - r.mean()) / r.std()
                rng.shuffle(E)

                p_cf, F_dual, s_star = closed_form(E, alpha, T)
                F_cf, _, H_cf = F_of_p(E, p_cf, alpha, T)
                p_md = projected_gd(E, alpha, T)
                F_md, _, H_md = F_of_p(E, p_md, alpha, T)

                t0 = time.time()
                p_q, cv_q, H_q, _c = run_cvar_vqe(E, alpha, T, n=n, layers=layers,
                                                  iters=80, restarts=1, seed=0)
                dt = time.time() - t0
                F_q = float(cv_q - T * H_q)

                # how well does the circuit's p approximate the closed-form optimum?
                tv = 0.5 * float(np.abs(p_q - p_cf).sum())
                rows.append(dict(n=n, alpha=alpha, T=T,
                                 F_closed=F_cf, F_dual=F_dual, F_mirror=F_md, F_vqe=F_q,
                                 gap_vqe=F_q - F_cf, dualgap=abs(F_cf - F_dual),
                                 H_closed=H_cf, H_vqe=H_q, TV=tv, s_star=s_star,
                                 vqe_sec=dt))

    hdr = ("  n alpha    T |    F_closed      F_dual     F_mirror        F_vqe |"
           "  VQE-CF gap   dualgap |  H_cf   H_vqe |   TV  | vqe s")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['n']:3d} {r['alpha']:5.2f} {r['T']:4.2f} |"
              f" {r['F_closed']:11.6f} {r['F_dual']:11.6f} {r['F_mirror']:12.6f}"
              f" {r['F_vqe']:12.6f} | {r['gap_vqe']:+11.6f} {r['dualgap']:9.2e} |"
              f" {r['H_closed']:6.3f} {r['H_vqe']:7.3f} | {r['TV']:.3f} | {r['vqe_sec']:5.1f}")

    print()
    print("KEY: F_closed  = hinged-Gibbs closed form (primal value)")
    print("     F_dual    = 1-D concave dual value  (equality => strong duality holds)")
    print("     F_mirror  = independent mirror-descent optimum over the full simplex")
    print("     F_vqe     = shipped run_cvar_vqe (80 Adam iters, exact param-shift grad)")
    print("     gap_vqe>0 => the circuit is STRICTLY WORSE than the closed form")


if __name__ == "__main__":
    main()
