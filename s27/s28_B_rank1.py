#!/usr/bin/env python
"""s27/s28_B_rank1.py -- mechanism check for F5: is the hopping term's gradient-variance decay
the rank-one (Perron) component of A?

hop-only cost -<psi|A|psi> measured as in `s28_B_train.py` with A replaced by
  (a) A itself (the F5 row, reproduced),
  (b) its rank-one part  v1 v1^T  (v1 the unit Perron vector; A has unit spectral norm),
  (c) the residual  A - v1 v1^T,
  (d) a diagonal observable with the SAME spectrum as A (the eigenvalues of A on the basis
      states, a random assignment): the control that separates "off-diagonal" from "spectrum".
Also the Perron vector's participation ratio and its overlap with the uniform state. Native-free,
no RMSD. Writes `s28_B_rank1.json`.
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
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402

OUT = os.path.join(B.RESULTS, "s28_B_rank1.json")
NS = [4, 5, 6, 7, 8, 9]
N_THETA = 120


def var_hop(n, A, seed=1009):
    circ = Q.StatevectorCircuit(n, B.LAYERS)
    rng = np.random.default_rng(seed)
    g0 = []
    for _ in range(N_THETA):
        th = rng.normal(0.0, 0.6, circ.n_params())
        g0.append(float(B.grad_hop_paramshift(circ, th, A)[0]))
    return float(np.var(g0, ddof=1))


def main():
    pdbs = P.targets()
    pick = pdbs[::11][:12]
    rows = []
    t0 = time.time()
    for pdb in pick:
        cand, ch, _ = RP.channels_for(pdb)
        E_full = RP.zr(ch["DIS"])
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        order = RP.topm(E_full, cand.k, key)
        D = B.pairwise_rmsd_matrix(cand.W)
        for n in NS:
            dim = 1 << n
            if dim <= cand.k:
                g = B.kernel_graph(D[np.ix_(order[:dim], order[:dim])])
                A = g["A"]
            else:
                A = B.pad_graph(B.kernel_graph(D)["A"], dim)
            w, v = np.linalg.eigh(A)
            v1 = v[:, -1] * np.sign(v[:, -1].sum())
            R1 = np.outer(v1, v1) * w[-1]
            res = A - R1
            rng = np.random.default_rng(7)
            diag_same_spectrum = np.diag(w[rng.permutation(dim)])
            u = np.full(dim, 1 / np.sqrt(dim))
            rows.append(dict(pdb=pdb, n=n, dim=dim, lam1=float(w[-1]), lam2=float(w[-2]),
                             perron_pr=float(1.0 / np.sum(v1 ** 4)), perron_uniform_overlap=float((u @ v1) ** 2),
                             var_A=var_hop(n, A), var_rank1=var_hop(n, R1), var_residual=var_hop(n, res),
                             var_diag_same_spectrum=var_hop(n, diag_same_spectrum)))
        print(f"  {pdb} ({(time.time()-t0)/60:.1f} min)", flush=True)
    summ = {}
    for k in ("var_A", "var_rank1", "var_residual", "var_diag_same_spectrum"):
        med = [float(np.median([r[k] for r in rows if r["n"] == n])) for n in NS]
        slope = float(np.polyfit(np.array(NS, float), np.log2(np.array(med)), 1)[0])
        summ[k] = dict(median_per_n={str(n): m for n, m in zip(NS, med)}, log2_slope_per_qubit=slope)
    for k in ("perron_pr", "perron_uniform_overlap", "lam2"):
        summ[k] = {str(n): float(np.median([r[k] for r in rows if r["n"] == n])) for n in NS}
    ST.save_atomic(OUT, dict(rows=rows, summary=summ, n_theta=N_THETA, targets=pick), module_file=__file__)
    print(f"\n{'observable':24s} " + " ".join(f"{'n='+str(n):>11}" for n in NS) + f" {'slope':>8}")
    for k in ("var_A", "var_rank1", "var_residual", "var_diag_same_spectrum"):
        s = summ[k]
        print(f"{k:24s} " + " ".join(f"{s['median_per_n'][str(n)]:11.3e}" for n in NS) + f" {s['log2_slope_per_qubit']:+8.3f}")
    print(f"{'perron PR / dim':24s} " + " ".join(f"{summ['perron_pr'][str(n)]/(1<<n):11.3f}" for n in NS))
    print(f"{'perron |<u|v1>|^2':24s} " + " ".join(f"{summ['perron_uniform_overlap'][str(n)]:11.3f}" for n in NS))
    print(f"{'lambda_2 / lambda_1':24s} " + " ".join(f"{summ['lam2'][str(n)]:11.3f}" for n in NS))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
