#!/usr/bin/env python
"""s27/s28_B2_rank1.py -- S28-L23 caveat (a): decompose the kNN graph's hop-only gradient
variance into its rank-one (Perron) part and the remainder, as `s28_B_rank1.py` did for the
Gaussian graph. Same draws, same targets, same registers. Native-free, no RMSD.
Writes `s28_B2_rank1.json`.
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

from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402
from s27 import s28_B_rank1 as R1          # noqa: E402
from s27 import s28_B2_knn as B2           # noqa: E402

OUT = os.path.join(B.RESULTS, "s28_B2_rank1.json")
ROWS = os.path.join(B.RESULTS, "s28_B2_rank1_rows.jsonl")
NS = [4, 5, 6, 7, 8, 9]


def main():
    pdbs = P.targets()
    pick = pdbs[::11][:12]
    done = B._done_pdbs(ROWS)
    t0 = time.time()
    for pdb in pick:
        if pdb in done:
            continue
        cand, ch, _ = RP.channels_for(pdb)
        E_full = RP.zr(ch["DIS"])
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        order = RP.topm(E_full, cand.k, key)
        D = B.pairwise_rmsd_matrix(cand.W)
        rows = []
        for k in B2.KS:
            for n in NS:
                dim = 1 << n
                if dim <= cand.k:
                    g = B2.knn_graph(D[np.ix_(order[:dim], order[:dim])], k)
                    A = g["A"]
                else:
                    g = B2.knn_graph(D, k)
                    A = B.pad_graph(g["A"], dim)
                w, v = np.linalg.eigh(A)
                v1 = v[:, -1] * np.sign(v[:, -1].sum() or 1.0)
                Rk = np.outer(v1, v1) * w[-1]
                res = A - Rk
                rows.append(dict(pdb=pdb, k=k, n=n, dim=dim, lam1=float(w[-1]), lam2=float(w[-2]),
                                 n_components=g["n_components"],
                                 perron_pr_over_dim=float(1.0 / np.sum(v1 ** 4) / dim),
                                 perron_uniform_overlap=float((np.full(dim, dim ** -0.5) @ v1) ** 2),
                                 var_A=R1.var_hop(n, A), var_rank1=R1.var_hop(n, Rk), var_residual=R1.var_hop(n, res)))
        with open(ROWS, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  {pdb} ({(time.time()-t0)/60:.1f} min)", flush=True)
    rows = B.load_rows(ROWS)
    summ = {}
    for k in B2.KS:
        for key in ("var_A", "var_rank1", "var_residual"):
            med = [float(np.median([r[key] for r in rows if r["k"] == k and r["n"] == n])) for n in NS]
            slope = float(np.polyfit(np.array(NS, float), np.log2(np.array(med)), 1)[0])
            summ[f"k{k}|{key}"] = dict(median_per_n={str(n): m for n, m in zip(NS, med)}, log2_slope_per_qubit=slope)
        summ[f"k{k}|residual_share"] = {str(n): float(np.median([r["var_residual"] / r["var_A"] for r in rows if r["k"] == k and r["n"] == n])) for n in NS}
        summ[f"k{k}|rank1_share"] = {str(n): float(np.median([r["var_rank1"] / r["var_A"] for r in rows if r["k"] == k and r["n"] == n])) for n in NS}
        summ[f"k{k}|n_components_max"] = {str(n): int(max(r["n_components"] for r in rows if r["k"] == k and r["n"] == n)) for n in NS}
    ST.save_atomic(OUT, dict(rows=rows, summary=summ, targets=pick), module_file=__file__)
    print(f"\n{'observable':26s} " + " ".join(f"{'n='+str(n):>11}" for n in NS) + f" {'slope':>8}")
    for k in B2.KS:
        for key in ("var_A", "var_rank1", "var_residual"):
            s = summ[f"k{k}|{key}"]
            print(f"{('k%d ' % k) + key:26s} " + " ".join(f"{s['median_per_n'][str(n)]:11.3e}" for n in NS) + f" {s['log2_slope_per_qubit']:+8.3f}")
        print(f"{('k%d ' % k) + 'residual share':26s} " + " ".join(f"{summ[f'k{k}|residual_share'][str(n)]:11.3f}" for n in NS))
        print(f"{('k%d ' % k) + 'rank1 share':26s} " + " ".join(f"{summ[f'k{k}|rank1_share'][str(n)]:11.3f}" for n in NS))
        print(f"{('k%d ' % k) + 'max components':26s} " + " ".join(f"{summ[f'k{k}|n_components_max'][str(n)]:11d}" for n in NS))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
