"""SPRINT 15 / QGEOM -- PART C(iv), CORRECTED: cost locality with the term count CONTROLLED.

`qgeom_grad.locality` builds an objective as a sum of `terms=32` random Pauli-Z products of
exact weight `w`.  At `w=1` there are only `n` possible supports and at `w=n` only ONE, so
those two columns are built from far fewer distinct terms than the middle of the range.  The
resulting profile is a deep U -- enormous at both ends, flat in the middle -- and the "slope
in w" fitted through it is meaningless.  **That is a confound in my own first measurement and
this module removes it.**

CORRECTED DESIGN
    * every weight uses the SAME number of distinct terms, and only weights with
      `C(n, w) >= terms` are compared;
    * the objective's variance is exactly 1 at every weight by construction (distinct
      Z-strings are orthogonal and +/-1 valued), so scale is not a confound;
    * the endpoints are still reported, labelled as term-count-limited, because their values
      are real and interesting -- they are just not part of the locality trend.

THE QUESTION.  Sprint 14 records that cost locality does NOT explain trainability at 6-18
qubits ("the measured ansatz kernel is flat in Pauli weight; what matters is n").  This tests
that on a controlled family, and reports the locality effect and the WIDTH effect in the same
units (log2 of gradient variance) so they can be compared directly.

    python -m s15.qgeom_loc
"""
from __future__ import annotations

from math import comb

import numpy as np

from s15 import qgeom_lib as G
from s15.qgeom_grad import grad_var

TAG = "loc"


def pauli_weight_objective(nq, w, terms, seed=0):
    """Sum of exactly `terms` DISTINCT random Pauli-Z products of weight `w`, variance 1."""
    r = np.random.default_rng(1000 * seed + 10 * nq + w)
    idx = np.arange(1 << nq, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(nq - 1, -1, -1)[None, :]) & 1).astype(np.int8)
    seen, out = set(), np.zeros(1 << nq)
    guard = 0
    while len(seen) < terms and guard < 100 * terms:
        guard += 1
        s = tuple(sorted(r.choice(nq, size=w, replace=False)))
        if s in seen:
            continue
        seen.add(s)
        out += (1 - 2 * (bits[:, list(s)].sum(1) % 2)).astype(float)
    return out / np.sqrt(len(seen)), len(seen)


def main(ns=(8, 10, 12, 14), terms=28, pat="ring", L=2, seeds=64):
    print("=" * 112)
    print("C4-CORRECTED. COST LOCALITY AT FIXED WIDTH, TERM COUNT CONTROLLED")
    print(f"    every weight uses exactly {terms} distinct Pauli-Z strings; variance = 1")
    print("    at every weight; only weights with C(n,w) >= terms enter the trend.")
    print("=" * 112)
    out = {}
    for n in ns:
        ws = list(range(1, n + 1))
        vs, ok = [], []
        for w in ws:
            E, got = pauli_weight_objective(n, w, terms)
            c = G.FlexCircuit(n, L, pat, 2)
            vs.append(grad_var(c, E, 1.0, seeds)["var_g1"])
            ok.append(got >= terms and comb(n, w) >= terms)
        trend = [(w, v) for w, v, o in zip(ws, vs, ok) if o]
        tw = np.array([t[0] for t in trend], float)
        tv = np.log2(np.maximum([t[1] for t in trend], 1e-300))
        sl = float(np.polyfit(tw, tv, 1)[0]) if len(tw) > 2 else float("nan")
        rng = float(tv.max() - tv.min()) if len(tv) else float("nan")
        out[f"n{n}"] = {"weights": ws, "var_g1": vs, "in_trend": ok,
                        "trend_weights": [int(x) for x in tw],
                        "log2_range_in_trend": rng, "slope_log2_per_weight": sl,
                        "log2_range_all": float(np.log2(max(vs) / max(min(vs), 1e-300)))}
        print(f"\n  n={n}   (weights in the trend: {[int(x) for x in tw]})")
        print("    " + "  ".join(f"w{w}{'' if o else '*'}={v:.2e}"
                                 for w, v, o in zip(ws, vs, ok)))
        print(f"    trend: log2 range {rng:.2f}  slope {sl:+.3f} per unit weight   "
              f"(all weights incl. endpoints: {out[f'n{n}']['log2_range_all']:.2f})")
        G.ck(TAG, "C4c_locality_controlled", out)
    print()
    print("  * = term-count-limited (C(n,w) < terms); reported, not in the trend.")
    print()
    print("  COMPARISON IN THE SAME UNITS.  C1 measures the WIDTH effect at ~ -0.23 log2")
    print("  per qubit, i.e. ~2.8 log2 over 6->18 qubits.  The locality effect over the")
    print("  well-populated weight range is printed above.  Whichever is larger is the")
    print("  axis that governs trainability at this scale.")
    return out


if __name__ == "__main__":
    main()
