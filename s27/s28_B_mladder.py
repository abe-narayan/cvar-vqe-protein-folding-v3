#!/usr/bin/env python
"""s27/s28_B_mladder.py -- the S27 T9 decomposition for every S28-B R1 row:

    RMSD_R1 = RMSD_top-75 + [RMSD_top-m - RMSD_top-75] + eps,   eps = RMSD_R1 - RMSD_top-m(E)

At J = 0, eps is the theorem reading back (<= 6.3e-4 on tie-free pools, S27 L6). At J > 0 the
tail is still a subset of an E-prefix (the operator's property), so eps is non-zero only through
holes; the m-ladder term is what a changed rung buys. Reads `s28_B_rows.jsonl`, writes
`s28_B_mladder.json`. ORACLE: reads natives for the RMSDs.
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402

OUT = os.path.join(B.RESULTS, "s28_B_mladder.json")


def main():
    rows = B.load_rows()
    by_pdb = {}
    for r in rows:
        by_pdb.setdefault(r["pdb"], []).append(r)
    out_rows = []
    for pdb, rs in sorted(by_pdb.items()):
        cand, ch, _ = RP.channels_for(pdb)
        E = RP.zr(ch["DIS"])
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        order = RP.topm(E, cand.k, key)
        cache = {}

        def topm_rmsd(m):
            if m not in cache:
                C, _ = I.coordinate_average(cand.W[order[:m]])
                cache[m] = float(I.ca_rmsd(C, cand.nat_ca))
            return cache[m]

        r75 = topm_rmsd(75)
        for r in rs:
            m = int(r["m"])
            rm = topm_rmsd(m) if m > 0 else float("nan")
            out_rows.append(dict(pdb=pdb, arm=B.arm_of(r), source=r["source"], seed=r["seed"],
                                 graph=r["graph"], J=r["J"], m=m, rmsd_R1=r["R1"]["rmsd"],
                                 rmsd_topm=rm, rmsd_top75=r75, eps=float(r["R1"]["rmsd"] - rm),
                                 mladder=float(rm - r75), n_holes=int(r["gate_n_holes"])))
    summ = {}
    arms = sorted({r["arm"] for r in out_rows})
    for a in arms:
        rs = [r for r in out_rows if r["arm"] == a]
        summ[a] = dict(n=len(rs), m_mean=float(np.mean([r["m"] for r in rs])),
                       eps_max_abs=float(np.max(np.abs([r["eps"] for r in rs]))),
                       eps_mean=float(np.mean([r["eps"] for r in rs])),
                       mladder_mean=float(np.mean([r["mladder"] for r in rs])),
                       holes_max=int(max(r["n_holes"] for r in rs)),
                       frac_with_holes=float(np.mean([r["n_holes"] > 0 for r in rs])))
    ST.save_atomic(OUT, dict(rows=out_rows, summary=summ), module_file=__file__)
    print(f"{'arm':22s} {'m':>6} {'eps max':>9} {'eps mean':>9} {'mladder':>8} {'holes':>6} {'f_holes':>7}")
    for a in arms:
        s = summ[a]
        print(f"{a:22s} {s['m_mean']:6.1f} {s['eps_max_abs']:9.2e} {s['eps_mean']:+9.4f} {s['mladder_mean']:+8.4f} {s['holes_max']:6d} {s['frac_with_holes']:7.3f}")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
