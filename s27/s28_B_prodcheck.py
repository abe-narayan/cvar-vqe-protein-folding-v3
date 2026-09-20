#!/usr/bin/env python
"""s27/s28_B_prodcheck.py -- the production comparator re-projected IN LANE B's OWN PROCESS on
all 126 targets, compared bit-exactly with S27's `chain_rows.jsonl :: DIS`.

Why (S28-L18 / S28-L27b / S28-L38): lane B's chain rows (`s28_B_hop.py :: chain_main`) and S27's
production chain rows (`run_vqe_chain.py --chain`) are two jobs; lane D reads a contrast across
jobs as "crossing a code path" and attaches the projection's branch-flip floor (mean 0.006 A,
tail 0.5 A on one target). Both jobs call `s24.d_harness.readout_projected` on the same kind of
cloud, and S28-L7 reproduced S27's rows bit-exactly on 6 targets. This script does it on 126:
the DIS top-75 uniform average (`readout_uniform`, the production readout) projected through
the SAME call lane B's chain rows use, then `rmsd_chain` compared with S27's stored value. If
the 126 values are identical to the last digit, every "vs production" contrast in lane B's
verdict is on one code path and the floor does not enter.

ORACLE: `rmsd_chain` reads `nat_ca` (post-hoc scoring only; the cloud never sees a native).
Writes `s27/results/s28_B_prodcheck_rows.jsonl` (per target, resumable) and
`s27/results/s28_B_prodcheck.json`.
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

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import run_vqe_chain as RV        # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402

ROWS = os.path.join(B.RESULTS, "s28_B_prodcheck_rows.jsonl")
OUT = os.path.join(B.RESULTS, "s28_B_prodcheck.json")
S27_CHAIN = os.path.join(B.RESULTS, "chain_rows.jsonl")


def s27_dis_rows():
    ref = {}
    for r in B.load_rows(S27_CHAIN):
        if r["config"] == "DIS":
            ref[r["pdb"]] = r
    return ref


def project_production(pdb: str):
    """S27's production chain row, recomputed: DIS top-75 (stable key) -> uniform average ->
    `readout_projected` (the call `s28_B_hop.chain_main` makes) -> CA-RMSD to the native."""
    cand, ch, _ = RP.channels_for(pdb)
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    E = RV.energy_for("DIS", ch, pdb, key)
    top = np.lexsort((key, E))[:RP.M]
    C, _ = H.readout_uniform(cand, top)
    ca = H.readout_projected(cand, C)
    return dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold),
                rmsd_cloud=float(I.ca_rmsd(C, cand.nat_ca)),
                rmsd_chain=float(I.ca_rmsd(ca, cand.nat_ca)))


def main():
    ref = s27_dis_rows()
    done = B._done_pdbs(ROWS)
    t0 = time.time()
    for i, pdb in enumerate(P.targets()):
        if pdb in done:
            continue
        t1 = time.time()
        r = project_production(pdb)
        r["ref_cloud"] = float(ref[pdb]["rmsd_cloud"])
        r["ref_chain"] = float(ref[pdb]["rmsd_chain"])
        r["d_cloud"] = abs(r["rmsd_cloud"] - r["ref_cloud"])
        r["d_chain"] = abs(r["rmsd_chain"] - r["ref_chain"])
        r["identical_chain"] = bool(r["rmsd_chain"] == r["ref_chain"])
        r["secs"] = float(time.time() - t1)
        with open(ROWS, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
        print(f"  [{i+1}/126] {pdb} chain {r['rmsd_chain']:.6f} ref {r['ref_chain']:.6f} "
              f"d {r['d_chain']:.2e} ({r['secs']:.1f}s, elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    rows = B.load_rows(ROWS)
    dc = np.array([r["d_chain"] for r in rows]); dk = np.array([r["d_cloud"] for r in rows])
    summ = dict(n=len(rows), n_identical_chain=int(sum(r["identical_chain"] for r in rows)),
                n_identical_cloud=int(sum(r["d_cloud"] == 0.0 for r in rows)),
                max_d_chain=float(dc.max()), mean_d_chain=float(dc.mean()),
                max_d_cloud=float(dk.max()),
                worst_pdb=rows[int(np.argmax(dc))]["pdb"],
                mean_chain_recomputed=float(np.mean([r["rmsd_chain"] for r in rows])),
                mean_chain_s27=float(np.mean([r["ref_chain"] for r in rows])))
    ST.save_atomic(OUT, dict(kind="production DIS top-75 chain rows re-projected in lane B's process vs S27 chain_rows DIS",
                             rows=rows, summary=summ), module_file=__file__)
    print(json.dumps(summ, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
