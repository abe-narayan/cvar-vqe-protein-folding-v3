#!/usr/bin/env python
"""s27/s28_D_anchor_B.py -- lane D: does lane B's J = 0 loop reproduce S27's DIS VQE row bit-for-bit?

For the named targets (default 1A13, 2N9M): run `s27.s28_B_hop.run_hop_vqe` at J = 0 with the
zero graph, read the tail with `readout_tail`, score the uniform average (ORACLE, post hoc),
and compare with (a) `s24.d_harness.arm_vqe` on the same energies (the p vector must be
identical, `==` on every float) and (b) the stored `s27/results/vqe_rows.jsonl` row
(config DIS, seed 0/1). Writes `s27/results/s28_D_anchor_B.json`.

Usage:  python s27/s28_D_anchor_B.py [--pdbs 1A13 2N9M]
"""
from __future__ import annotations

import argparse
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

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbs", nargs="*", default=["1A13", "2N9M"])
    a = ap.parse_args()
    stored = {}
    with open(os.path.join(HERE, "results", "vqe_rows.jsonl"), encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            if r["config"] == "DIS":
                stored[(r["pdb"], int(r["seed"]))] = r
    out = {"targets": {}}
    for pdb in a.pdbs:
        cand, ch, _ = RP.channels_for(pdb)
        E = RP.zr(ch["DIS"])
        res = {}
        for seed in (0, 1):
            t0 = time.time()
            enc = B.QC.Encoding(E)
            zero = np.zeros((enc.dim, enc.dim))
            p_b, cv_b, H_b, hv_b, F_b, th_b, circ = B.run_hop_vqe(enc.E, zero, 0.0, seed=seed, n=enc.n_qubits)
            r1 = B.readout_tail(enc, p_b, B.ALPHA)
            rmsd_b = P.rmsd_of_set(cand, r1["idx"])
            v = H.arm_vqe(cand, E, alpha=0.18, T=0.5, layers=3, iters=80, seed=seed)
            # the harness's p is not returned; recompute it the harness's way
            p_h, cv_h, H_h, _ = B.Q.run_cvar_vqe(enc.E, 0.18, 0.5, n=enc.n_qubits, layers=3, iters=80,
                                                 restarts=1, seed=seed, lr=0.15)
            rmsd_h = P.rmsd_of_set(cand, v["cands"])
            st = stored.get((pdb, seed))
            res[str(seed)] = dict(
                p_identical=bool(np.array_equal(p_b, p_h)),
                p_max_abs_diff=float(np.abs(p_b - p_h).max()),
                cvar_b=float(cv_b), cvar_h=float(cv_h),
                m_b=int(r1["m"]), m_h=int(v["m"]),
                idx_identical=bool(np.array_equal(np.sort(r1["idx"]), np.sort(v["cands"]))),
                rmsd_b=float(rmsd_b), rmsd_h=float(rmsd_h),
                rmsd_stored=float(st["rmsd_vqe"]) if st else None,
                m_stored=int(st["m"]) if st else None,
                abs_diff_vs_stored=abs(float(rmsd_b) - float(st["rmsd_vqe"])) if st else None,
                secs=time.time() - t0)
            print(pdb, "seed", seed, json.dumps(res[str(seed)]))
        out["targets"][pdb] = res
    out["all_identical"] = bool(all(r["p_identical"] and r["idx_identical"] and (r["abs_diff_vs_stored"] == 0.0)
                                    for t in out["targets"].values() for r in t.values()))
    print("ALL IDENTICAL:", out["all_identical"])
    ST.save_atomic(os.path.join(HERE, "results", "s28_D_anchor_B.json"), out,
                   complete_keys=["targets", "all_identical"], module_file=__file__)


if __name__ == "__main__":
    main()
