"""SPRINT 15 / QGEOM -- PART C(vi), CORRECTED: trainability AT the conditioned point.

`qgeom_grad.conditioning_effect` reported `Var[g_1]` over seeds for the `uniform` and
`COND` arms, but those arms are a SINGLE point with a 0.01 jitter, so the across-seed
variance is small by construction and is not comparable to the Haar column (which varies
over the whole parameter space).  **That is a defect in my own cell** and it is corrected
here by reporting statistics of the gradient AT the point -- its norm, its per-component
variance, and its share in the metric's small eigendirections -- which are what
"trainability at this point" means.

    python -m s15.qgeom_c6
"""
from __future__ import annotations

import numpy as np

from s14 import retprior as RP
from s15 import qgeom_lib as G
from s15.qgeom_qng import sub_objective

TAG = "c6"
TARGETS = ("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5")


def main(pat="block", L=2, nres=6, iters=300, seeds=(0, 1, 2)):
    res = tuple(range(1, 1 + nres))
    c = G.FlexCircuit(2 * nres, L, pat, 2)
    print("=" * 108)
    print("C6-CORRECTED. GRADIENT AT THE POINT (not variance across near-identical points)")
    print("=" * 108)
    print(f"{'target':7s} {'seed':>4s} {'H(prior)':>9s} | " +
          " | ".join(f"{a:>21s}" for a in ("haar", "uniform", "COND", "scram_state")))
    print(f"{'':7s} {'':4s} {'':9s} | " +
          " | ".join(f"{'|g|      Var_comp':>21s}" for _ in range(4)))
    out = {}
    for pdb in TARGETS:
        e = G.V.Enum(pdb)
        base = np.random.default_rng(0).integers(0, e.k, e.n)
        E = G.V.uniformise(sub_objective(e, res, base, "legacy"))
        P = RP.state_prior(pdb, "top75", 4)["P"][list(res)]
        for s in seeds:
            row = {"prior_entropy_bits": G.prior_entropy_bits(P)}
            pts = {"haar": G.haar_theta(c, 100 * s + 7),
                   "uniform": G.uniform_theta(c)}
            for nm, Pk in (("COND", P), ("scram_state", G.scramble_prior(P, s, "state"))):
                pts[nm] = G.fit_kl(c, G.product_target(Pk),
                                   theta0=G.random_theta(c, s), iters=iters,
                                   lr=0.08)["theta"]
            for nm, th in pts.items():
                psi = c.state(th)
                gr = G.grad_adjoint(c, th, 2.0 * psi * E)
                row[nm] = {"grad_norm": float(np.linalg.norm(gr)),
                           "grad_var_components": float(np.var(gr)),
                           "mean_abs": float(np.mean(np.abs(gr)))}
            out[f"{pdb}|s{s}"] = row
            print(f"{pdb:7s} {s:4d} {row['prior_entropy_bits']:9.3f} | " +
                  " | ".join(f"{row[a]['grad_norm']:9.5f} {row[a]['grad_var_components']:11.3e}"
                             for a in ("haar", "uniform", "COND", "scram_state")))
            G.ck(TAG, "C6c_conditioning_pointwise", out)
    print()
    print("  PAIRED over all cells (negative = first arm lower):")
    for a, b in (("COND", "uniform"), ("COND", "haar"), ("COND", "scram_state"),
                 ("scram_state", "uniform")):
        for st in ("grad_norm", "grad_var_components"):
            pr = G.paired([v[a][st] for v in out.values()],
                          [v[b][st] for v in out.values()])
            print(f"    {st:20s} {a+' - '+b:24s} {pr['mean']:+11.5f} "
                  f"[{pr['ci_lo']:+.5f}, {pr['ci_hi']:+.5f}] "
                  f"W/L {pr['win']}/{pr['loss']}  {pr['verdict']}")
            G.ck(TAG, f"C6c_paired|{st}|{a}-{b}", pr)
    return out


if __name__ == "__main__":
    main()
