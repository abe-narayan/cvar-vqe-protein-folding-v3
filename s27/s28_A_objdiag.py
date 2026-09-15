#!/usr/bin/env python
"""s27/s28_A_objdiag.py -- ORACLE DIAGNOSTIC of the native-free objective term (PREREG section 6).

For every target: the surrogate S~ and the shipped S of (a) the NATIVE posed in the frame
(ORACLE), (b) the ORACLE circuit optimum, (c) production's average, (d) every pool member, so
the native's percentile under the objective inside the pool is known.  Answers "does the
objective even prefer the near-native structure?" before any recognition arm is read.
Writes `s27/results/s28_A_objdiag.json`.  Everything here is ORACLE except the S values of
production and the pool members, which are native-free quantities reported beside them.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A_amp as A             # noqa: E402


def main():
    from s25 import phys_lib as P
    rows = []
    for k, pdb in enumerate(P.targets()):
        cand, dis, top, dg = A.load_pool(pdb)
        frame = A.Frame(cand.W, top)
        sur = A.Surrogate(dg, cand.n)
        z = np.load(os.path.join(A.STRUCTS, f"{pdb}.npz"))
        natp = I.superpose_batch(np.asarray(cand.nat_ca, float)[None], frame.ref)[0]   # ORACLE
        s_nat = float(sur.value_grad(natp)[0])
        s_orc = float(sur.value_grad(np.asarray(z["oracle_circ"], float))[0])
        s_prod = float(sur.value_grad(np.asarray(z["prod"], float))[0])
        s_pool = np.array([sur.value_grad(w)[0] for w in frame.Wp])
        rr = np.asarray(cand.oracle_rr, float)
        rows.append(dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold),
                         S_native=s_nat, S_oracle_circ=s_orc, S_prod=s_prod,
                         S_pool_min=float(s_pool.min()), S_pool_mean=float(s_pool.mean()),
                         native_pctile_in_pool=float((s_pool < s_nat).mean()),
                         S_shipped_native=sur.shipped(natp), S_shipped_prod=sur.shipped(np.asarray(z["prod"], float)),
                         S_at_pool_best=float(s_pool[int(np.argmin(rr))]),
                         prod_below_native=bool(s_prod < s_nat), prod_below_pool_min=bool(s_prod < s_pool.min())))
        if (k + 1) % 20 == 0:
            print(f"  {k+1}/126", flush=True)
    S = lambda key: np.array([r[key] for r in rows], float)
    out = dict(rows=rows,
               mean_S_native=float(S("S_native").mean()), mean_S_oracle_circ=float(S("S_oracle_circ").mean()),
               mean_S_prod=float(S("S_prod").mean()), mean_S_pool_min=float(S("S_pool_min").mean()),
               mean_S_pool_mean=float(S("S_pool_mean").mean()),
               native_pctile_in_pool_mean=float(S("native_pctile_in_pool").mean()),
               n_prod_below_native=int(sum(r["prod_below_native"] for r in rows)),
               n_prod_below_pool_min=int(sum(r["prod_below_pool_min"] for r in rows)),
               n_native_below_prod=int((S("S_native") < S("S_prod")).sum()))
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    ST.save_atomic(os.path.join(A.RESULTS, "s28_A_objdiag.json"), out, rows=rows, n_expected=126,
                   complete_keys=("pdb", "S_native", "S_prod", "S_pool_min"), module_file=__file__)


if __name__ == "__main__":
    main()
