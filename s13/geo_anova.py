"""SPRINT 13 GEO -- EXPERIMENT 2: how LOCAL is each energy model in torsion variables?

The BRIEF's section 4 is a claim about the Hamiltonian, and it is testable exactly:

    "In torsion space the chain builds sequentially, so a distance between residues i and j
     depends on every torsion between them.  A pairwise distance term is therefore NOT
     2-local in the torsion variables -- it is |i-j|-local. ... Legacy's terms are pairwise
     in CA/CB distance -> high-weight in torsion space; AMBER's are pairwise over ~230
     atoms -> higher-weight still, and denser."

Because the register is enumerated, the **exact functional (Sobol/ANOVA) decomposition** of
each energy model over the L per-residue torsion variables is available with no sampling
and no surrogate:

    E(s) = f0 + sum_i f_i(s_i) + sum_{i<j} f_ij(s_i,s_j) + sum_{i<j<l} f_ijl + ...

with the terms orthogonal under the uniform product measure, so
``Var(E) = V1 + V2 + V3 + V_{>=4}`` EXACTLY (the residual is float error, reported).

This is the measurement that decides whether "AMBER is harder" is a statement about
*interaction order* or only about *scale*.  Three quantities per model:

    * the order profile V1/V2/V3/V>=4 as shares of Var(E);
    * the effective interaction order -- smallest d with (V1+..+Vd)/Var >= 0.95;
    * the pair-coupling profile V2_ij against sequence separation |i-j|, which is the
      direct test of "|i-j|-local" against "2-local".

CONDITIONING -- CORRECTED 2026-09-05.  The headline is now the RANK-PRESERVING soft
compression (`geo_common.soft_compress`) applied identically to BOTH models, with raw and
99th-percentile-winsorised decompositions reported beside it and a `spike` diagnostic on
each.  The earlier headline used the 99th-percentile winsorisation alone, and that is NOT
enough: measured here, a winsorised AMBER table can still carry 0.99 of its variance in ten
configurations (1A13 L=4 k=4), and so can a raw LEGACY table (0.88 on 1A13 L=5 k=4).  A
decomposition of a table dominated by a handful of steric spikes reports the spikes, not the
force field.  See `s13/walsh_FINDINGS.md` (Pauli-spectrum agent) and section 4 of
`s13/geo_FINDINGS.md` for the correction record.

    python -m s13.geo_anova
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402

CELLS = [("1A13", 6, 4), ("1A1P", 6, 4), ("2BFI", 6, 4), ("1CEK", 6, 4),
         ("1DEP", 6, 4), ("1ID6", 6, 4),
         ("1A13", 7, 4), ("1A1P", 7, 4), ("2BFI", 7, 4),
         ("1A13", 12, 2), ("1A1P", 12, 2), ("2BFI", 12, 2),
         ("1A13", 4, 8), ("1A1P", 4, 8), ("2BFI", 4, 8)]


def _marg(T: np.ndarray, U):
    """Mean of T over every axis NOT in U, keeping U's axes in order."""
    L = T.ndim
    drop = tuple(a for a in range(L) if a not in U)
    return T.mean(axis=drop) if drop else T


def anova(E: np.ndarray, L: int, k: int, max_order: int = 3):
    """Exact Sobol decomposition up to `max_order`; returns variance shares."""
    T = np.asarray(E, float).reshape((k,) * L)
    V = float(T.var())
    f0 = float(T.mean())
    marg = {(): np.array(f0)}
    eff = {}
    Vd = {}
    pair = {}
    for d in range(1, max_order + 1):
        tot = 0.0
        for U in itertools.combinations(range(L), d):
            m = _marg(T, U)
            marg[U] = m
            f = m.copy()
            for r in range(d):
                for Vsub in itertools.combinations(U, r):
                    # broadcast the lower-order term back onto U's axes
                    sh = [k if a in Vsub else 1 for a in U]
                    src = marg[Vsub] if Vsub else np.array(f0)
                    f = f - (np.asarray(eff[Vsub]).reshape(sh) if Vsub
                             else np.asarray(src).reshape([1] * d))
            eff[U] = f
            v = float((f ** 2).mean())
            tot += v
            if d == 2:
                pair[U] = v
        Vd[d] = tot
    used = sum(Vd.values())
    return {
        "var_total": V,
        "share_order": {str(d): (Vd[d] / V if V > 0 else float("nan"))
                        for d in Vd},
        "share_ge_%d" % (max_order + 1): (V - used) / V if V > 0 else float("nan"),
        "cum_share": {str(d): (sum(Vd[e] for e in range(1, d + 1)) / V if V > 0
                               else float("nan")) for d in Vd},
        "pair_by_sep": _by_sep(pair, V, L),
    }


def _by_sep(pair, V, L):
    out = {}
    for (i, j), v in pair.items():
        out.setdefault(j - i, []).append(v)
    return {str(s): {"n": len(vs), "share_sum": float(np.sum(vs) / V),
                     "share_mean": float(np.mean(vs) / V)}
            for s, vs in sorted(out.items())}


def eff_order(a, thresh=0.95):
    for d in ("1", "2", "3"):
        if a["cum_share"].get(d, 0.0) >= thresh:
            return int(d)
    return 4


def main():
    rows = []
    for pdb, L, k in CELLS:
        p = os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")
        if not os.path.exists(p):
            print(f"  skip {pdb} L={L} k={k} (amber table not built yet)")
            continue
        r = {"pdb": pdb, "L": L, "k": k, "n_qubits": int(round(np.log2(k))) * L}
        for m in G.VARIANTS:
            E = G.variant_table(pdb, L, k, m)
            r[m] = {"soft": anova(G.soft_compress(E), L, k),
                    "winsor99": anova(G.winsor_hi(E, 0.99), L, k),
                    "raw": anova(E, L, k),
                    "spike_raw": G.spike_share(E),
                    "spike_winsor99": G.spike_share(G.winsor_hi(E, 0.99)),
                    "spike_soft": G.spike_share(G.soft_compress(E))}
            r[m]["eff_order_soft"] = eff_order(r[m]["soft"])
            r[m]["eff_order_winsor"] = eff_order(r[m]["winsor99"])
            r[m]["eff_order_raw"] = eff_order(r[m]["raw"])
        rows.append(r)
        print("  %-5s L=%-2d k=%d  SOFT 1/2/3/>=4  legacy %.3f %.3f %.3f %.3f | "
              "amber %.3f %.3f %.3f %.3f  (raw top10 spike leg %.3f amb %.3f)"
              % (pdb, L, k,
                 *[r["legacy_soft"]["soft"]["share_order"][d] for d in ("1", "2", "3")],
                 r["legacy_soft"]["soft"]["share_ge_4"],
                 *[r["amber_soft"]["soft"]["share_order"][d] for d in ("1", "2", "3")],
                 r["amber_soft"]["soft"]["share_ge_4"],
                 r["legacy"]["spike_raw"]["top10"], r["amber"]["spike_raw"]["top10"]),
              flush=True)
    G.write("geo_anova", {"what": "exact Sobol/ANOVA locality of the two energy models "
                                 "over the per-residue torsion variables",
                          "cells": rows})


if __name__ == "__main__":
    main()
