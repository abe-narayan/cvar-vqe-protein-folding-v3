"""SPRINT 15 / QGEOM -- PART A(iii): WHY QNG HELPS ON ONE COST AND NOT ANOTHER.

A2 measured the gradient's spectral profile in the eigenbasis of the metric FOR THE
EXPECTATION VALUE of a structural objective, found it depletes the small eigendirections
25x-60x below uniform, and I pre-registered from that: **QNG will be neutral at best**.

A3 then found QNG clearly HELPING on the KL distribution-fitting task in exactly the
ill-conditioned cells.  **My prediction is wrong on that cost**, and the obvious candidate
explanation is that the spectral profile is a property of the COST, not of the ansatz:

    grad_j <E>   = 2 sum_x E(x) psi(x) (d_j psi)(x)          -- weights states by amplitude
    grad_j KL    = -2 sum_x [t(x)/p(x)] psi(x) (d_j psi)(x)  -- weights by t/p, which BLOWS
                                                                UP on low-probability states

Low-probability states are exactly where the metric's small eigendirections live, so the KL
gradient should load the bottom of the spectrum and the expectation gradient should not.
This module measures both on identical circuits and identical points, so the explanation is
tested rather than asserted.

    python -m s15.qgeom_align
"""
from __future__ import annotations

import numpy as np

from core.quantum import cvar_exact
from s15 import qgeom_lib as G
from s15.qgeom_qng import LADDER, pack, sub_objective

TAG = "align"


def profile(circ, th, grad):
    g = G.fs_metric(circ, th)
    w, V = np.linalg.eigh((g + g.T) / 2.0)
    w = np.clip(w, 0, None)
    comp = (V.T @ grad) ** 2
    tot = max(comp.sum(), 1e-300)
    tol = 1e-9 * max(w.max(), 1e-30)
    pos = np.where(w > tol)[0]
    order = pos[np.argsort(w[pos])]
    nd = max(1, len(order) // 10)
    lam = w[order]
    share = comp[order] / tot
    # effective eigenvalue the gradient sees: the share-weighted harmonic-type mean
    lam_eff = float(1.0 / max((share / np.maximum(lam, 1e-30)).sum(), 1e-30))
    return {"bottom_decile": float(share[:nd].sum()),
            "top_decile": float(share[-nd:].sum()),
            "null_share": float(comp[w <= tol].sum() / tot),
            "cond": float(lam.max() / lam.min()) if len(lam) else np.nan,
            "lam_eff": lam_eff, "lam_max": float(lam.max()),
            "amplification": float(lam.max() / max(lam_eff, 1e-30))}


def main(n=12, seeds=6):
    print("=" * 116)
    print("A5. THE GRADIENT'S SPECTRAL PROFILE IS A PROPERTY OF THE COST, NOT THE ANSATZ")
    print("    Same circuits, same points; three costs.  `bottom decile` is the share of")
    print("    the gradient's squared norm in the metric's smallest 10% of eigenvalues")
    print("    (uniform would be 0.10).  QNG divides by the eigenvalue, so a gradient that")
    print("    LOADS the bottom decile is a gradient QNG can rescue -- and a gradient that")
    print("    AVOIDS it is one QNG can only amplify noise in.")
    print("=" * 116)
    pk = pack("1CS9")
    E = pk["E"]
    e = G.V.Enum("1CS9")
    tgt = G.V.uniformise(sub_objective(e, (1, 2, 3, 4, 5, 6), np.zeros(e.n, int), "rmsd"))
    tgt = np.exp(-6.0 * tgt)
    tgt = tgt / tgt.sum()
    costs = ("expectation", "cvar0.25", "KL")
    print(f"{'pattern':11s} {'L':>2s} {'cond':>9s} | " +
          " | ".join(f"{c:>22s}" for c in costs))
    print(f"{'':11s} {'':2s} {'':9s} | " +
          " | ".join(f"{'bot10  top10  amplif':>22s}" for c in costs))
    out = {}
    for pat, L in LADDER:
        c = G.FlexCircuit(n, L, pat, 2)
        rows = {k: [] for k in costs}
        conds = []
        for s in range(seeds):
            th = G.random_theta(c, 100 * s + 7)
            psi = c.state(th)
            p = psi ** 2
            p = p / p.sum()
            gs = {"expectation": 2.0 * (G.deriv_states(c, th) @ (psi * E))}
            gs["cvar0.25"] = G.grad_adjoint(
                c, th, 2.0 * psi * cvar_exact(E, p, 0.25)[2])
            gs["KL"] = -G.grad_adjoint(
                c, th, 2.0 * psi * (tgt / np.maximum(p, 1e-300)))
            for k in costs:
                pr = profile(c, th, gs[k])
                rows[k].append(pr)
                if k == "expectation":
                    conds.append(pr["cond"])
        out[f"{pat}|L{L}"] = {k: {kk: float(np.mean([r[kk] for r in rows[k]]))
                                  for kk in rows[k][0]} for k in costs}
        out[f"{pat}|L{L}"]["cond_median"] = float(np.median(conds))
        cells = []
        for k in costs:
            m = out[f"{pat}|L{L}"][k]
            cells.append(f"{m['bottom_decile']:6.4f} {m['top_decile']:6.3f} "
                         f"{m['amplification']:7.1f}")
        print(f"{pat:11s} {L:2d} {np.median(conds):9.2f} | " + " | ".join(cells))
        G.ck(TAG, "A5_cost_profiles", out)
    print()
    print("  `amplif` = lam_max / lam_eff, where lam_eff is the share-weighted eigenvalue")
    print("  the natural gradient actually divides by.  It is the factor by which QNG")
    print("  rescales this gradient relative to the stiffest direction -- large means QNG")
    print("  has real work to do, ~1 means the gradient already lives on the stiff")
    print("  directions and preconditioning can only move noise.")
    return out


if __name__ == "__main__":
    main()
