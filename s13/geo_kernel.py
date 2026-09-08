"""SPRINT 13 GEO -- EXPERIMENT 9: the ansatz gradient kernel v(w; n, depth), on its own.

WHY THIS EXPERIMENT EXISTS SEPARATELY
=====================================
`geo_pauli` measures the gradient-variance kernel

    v(w) = mean over |S| = w and over i of Var_theta[ d<Z_S>/dtheta_i ],   theta ~ U(-pi,pi)

as a by-product of the energy-model analysis, so it inherits that analysis's register sizes.
But **v contains no energy model at all** -- it is a property of the RY/CNOT ansatz and
nothing else, needs no AMBER call and no enumerated table, and can therefore be pushed to
whatever register size the simulator reaches.  That is exactly what the two open questions
need:

  Q1  Is the n-dependence of the gradient variance EXPONENTIAL (a barren plateau in the
      established sense) or POLYNOMIAL, over the range actually reachable?  Both fits, both
      R^2, at fixed Pauli weight and at fixed relative weight w/n.
  Q2  Does the OBSERVABLE'S LOCALITY control its gradient variance?  `geo_pauli` found
      v(w) non-monotone and nearly flat at n <= 12, contradicting the naive reading of the
      cost-locality framework.  Here the same measurement runs to n = 18-20 so the claim can
      be given an explicit regime of validity instead of an implicit one.

METHOD.  For each (n, layers) and each of `n_theta` random theta, one batched simulation
gives the exact probability Jacobian ``J = dp(x)/dtheta_i``; ONE Walsh-Hadamard transform of
its rows gives ``d<Z_S>/dtheta_i`` for ALL 2^n Pauli-Z strings simultaneously.  So the whole
Pauli basis is covered exactly, with no sampling over observables.

Exact throughout; no energy model, no AMBER, no native information.

    python -m s13.geo_kernel [max_n]
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from s13.geo_pauli import fwht, popcount              # noqa: E402

LAYERS = (1, 2, 3, 5, 8)
NS = (6, 8, 10, 12, 14, 16, 18)


def n_theta_for(n: int) -> int:
    return {6: 256, 8: 192, 10: 128, 12: 96, 14: 48, 16: 24, 18: 10, 20: 6}.get(n, 6)


def kernel(n: int, layers: int, n_theta: int | None = None, seed: int = 0):
    circ = G.circuit(n, layers=layers)
    P = circ.n_params()
    N = 1 << n
    nt = n_theta or n_theta_for(n)
    w = popcount(n)
    rng = np.random.default_rng(1_000_003 * seed + 1013 * n + layers)
    s1 = np.zeros((P, N)); s2 = np.zeros((P, N))
    for _ in range(nt):
        th = rng.uniform(-np.pi, np.pi, P)
        PR = circ.probs_batch(circ._shift_grid(th, np.pi / 2))
        Jh = fwht((PR[0::2] - PR[1::2]) / 2.0)
        s1 += Jh; s2 += Jh ** 2
    Var = s2 / nt - (s1 / nt) ** 2
    m = Var.mean(0)                                   # averaged over parameters
    vw = np.array([float(m[w == ww].mean()) if (w == ww).any() else 0.0
                   for ww in range(n + 1)])
    # per-layer profile: does depth position matter more than weight?
    byp = Var.mean(1).reshape(layers, n).mean(1)
    return {"n_qubits": n, "layers": layers, "P": P, "n_theta": nt,
            "v_of_weight": [float(x) for x in vw],
            "v_w1": float(vw[1]), "v_w2": float(vw[2]),
            "v_whalf": float(vw[max(1, n // 2)]), "v_wn": float(vw[n]),
            "v_max": float(vw[1:].max()), "argmax_weight": int(vw[1:].argmax() + 1),
            "v_wn_over_v_w1": float(vw[n] / max(vw[1], 1e-300)),
            "weight_decay_fit": G.loglog_fit(np.arange(1, n + 1), vw[1:]),
            "var_by_layer": [float(x) for x in byp]}


def main(max_n: int | None = None):
    ns = [n for n in NS if max_n is None or n <= int(max_n)]
    rows = []
    p = os.path.join(G.RESULTS, "geo_kernel.json")
    if os.path.exists(p):
        import json
        rows = json.load(open(p)).get("rows", [])
    done = {(r["n_qubits"], r["layers"]) for r in rows}
    for n in ns:
        for lay in LAYERS:
            if (n, lay) in done:
                continue
            t0 = time.time()
            r = kernel(n, lay)
            r["wall_s"] = time.time() - t0
            rows.append(r)
            print("  n=%-2d lay=%-2d P=%-3d nt=%-3d | v(1)=%.3e v(2)=%.3e v(n/2)=%.3e "
                  "v(n)=%.3e | v(n)/v(1)=%.3f argmax w=%d | %.0fs"
                  % (n, lay, r["P"], r["n_theta"], r["v_w1"], r["v_w2"], r["v_whalf"],
                     r["v_wn"], r["v_wn_over_v_w1"], r["argmax_weight"], r["wall_s"]),
                  flush=True)
            G.write("geo_kernel", {"what": "ansatz-only gradient-variance kernel over the "
                                          "complete Pauli-Z basis", "rows": rows,
                                   "scaling": scaling(rows)})
    G.write("geo_kernel", {"what": "ansatz-only gradient-variance kernel over the complete "
                                  "Pauli-Z basis", "rows": rows, "scaling": scaling(rows)})


def scaling(rows):
    """n-scaling at fixed weight and at fixed relative weight: exponential vs polynomial."""
    out = {}
    for lay in sorted({r["layers"] for r in rows}):
        sel = sorted([r for r in rows if r["layers"] == lay], key=lambda r: r["n_qubits"])
        if len(sel) < 3:
            continue
        ns = np.array([r["n_qubits"] for r in sel], float)
        for key in ("v_w1", "v_w2", "v_whalf", "v_wn"):
            v = np.array([r[key] for r in sel], float)
            out[f"lay{lay}_{key}"] = {"n": [int(x) for x in ns],
                                      "v": [float(x) for x in v], **G.loglog_fit(ns, v)}
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
