#!/usr/bin/env python
"""s27/run_trainability.py -- H8: does the Hamiltonian change the selector's trainability?

Gradient variance over theta ~ N(0, 0.6^2) of the deployed objective F = CVaR_alpha(E; p)
- T H(p) at the S25 suite's settings (n = 9 qubits over the 500 candidates + 12 padding
states at max + 10 sd, depth 3, alpha 0.18, T 0.5), exact parameter-shift gradient, for each
channel's energy vector under two monotone standardisations:

  zrank      the deployed currency: every channel becomes the same standardised rank ladder
             up to ties, so the variance is channel-independent BY CONSTRUCTION (measured to
             show it)
  asinh-MAD  a bounded standardisation that keeps the channel's own gaps

Property measurement on 12 targets (two per length class); no native is read.  Output:
`s27/results/trainability.json`.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q              # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import ham_lib as HL              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

CHANNELS = ["DIS", "LEG", "AMB", "CONS", "DSSPHB", "CONTACT", "DISTPOT", "RAMA", "CAGEO", "ENV"]
N_THETA, ALPHA, TEMP, LAYERS = 120, 0.18, 0.5, 3
OUT = os.path.join(RP.RESULTS, "trainability.json")


def padded(E, dim):
    E = np.asarray(E, float)
    pad = np.full(dim - len(E), E.max() + 10.0 * (E.std() + 1e-12))
    return np.concatenate([E, pad])


def main():
    pdbs = P.targets()
    pick = pdbs[::11][:12]
    circ = Q.StatevectorCircuit(9, LAYERS)
    Pn = circ.n_params()
    rows = []
    t0 = time.time()
    for pdb in pick:
        cand, ch, _ = RP.channels_for(pdb)
        for c in CHANNELS:
            for std_name, std in (("zrank", HL.zrank), ("asinh", HL.asinh_std)):
                E = padded(std(ch[c]), circ.dim)
                rng = np.random.default_rng(1009)                # same draws for every cell
                g0, gn = [], []
                for _ in range(N_THETA):
                    th = rng.normal(0.0, 0.6, Pn)
                    _, g, _, _, _ = Q.free_energy(circ, th, E, ALPHA, TEMP)
                    g0.append(float(g[0])); gn.append(float(np.dot(g, g)))
                g0 = np.asarray(g0); gn = np.asarray(gn)
                rows.append(dict(pdb=pdb, channel=c, std=std_name, var_g0=float(g0.var(ddof=1)),
                                 mean_sq_per_param=float(gn.mean() / Pn), E_sd=float(np.std(E[:cand.k])),
                                 E_range=float(np.ptp(E[:cand.k]))))
        print(f"  {pdb} done ({(time.time()-t0)/60:.1f} min)", flush=True)
    # summary per (channel, std): median over targets
    summ = {}
    for c in CHANNELS:
        for s in ("zrank", "asinh"):
            v = [r["var_g0"] for r in rows if r["channel"] == c and r["std"] == s]
            m = [r["mean_sq_per_param"] for r in rows if r["channel"] == c and r["std"] == s]
            summ[f"{c}|{s}"] = dict(var_g0_median=float(np.median(v)), var_g0_min=float(min(v)), var_g0_max=float(max(v)),
                                    mean_sq_per_param_median=float(np.median(m)))
    ref = summ["DIS|asinh"]["var_g0_median"]
    for k in summ:
        summ[k]["ratio_to_DIS_same_std"] = float(summ[k]["var_g0_median"] / summ[f"DIS|{k.split('|')[1]}"]["var_g0_median"])
    out = dict(targets=pick, n_theta=N_THETA, alpha=ALPHA, T=TEMP, n_qubits=9, layers=LAYERS, rows=rows, summary=summ)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"{'channel':10s} {'zrank var':>12} {'asinh var':>12} {'asinh/DIS':>10}")
    for c in CHANNELS:
        print(f"{c:10s} {summ[c+'|zrank']['var_g0_median']:12.4e} {summ[c+'|asinh']['var_g0_median']:12.4e} {summ[c+'|asinh']['ratio_to_DIS_same_std']:10.2f}")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
