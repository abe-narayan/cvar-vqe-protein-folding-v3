#!/usr/bin/env python
"""s27/s28_A2_prodcheck.py -- lane D's S28-L23 caveat (d), answered by measurement.

The A2 chain job (`s28A2_chain`, `_r2`) projected the 13 A2 arms; the production comparator on
the built chain is the PRIMARY chain job's projection (`s28A_chain_primary_1`) of the identical
stored cloud (`s27/results/s28_A_structs/<pdb>.npz :: prod`) through the identical code path
(`s28_A_amp.run_chain_target` -> `s12.instrument.project`).  The S28-L18 / S28-L27b floor needs an
INPUT difference (a branch flip under a 1e-13 perturbation); with the identical array it should
be bit-identical.  This re-projects `prod` on a few targets (2LNG is the 0.5 A branch-flip target
of S28-L27b) and compares to the stored rows.  Writes `s27/results/s28_A2_prodcheck.json`.
"""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A_amp as A             # noqa: E402

PDBS = ("2LNG", "1A13", "9BFL", "6RRO")


def main():
    rows = A._jsonl_done(os.path.join(A.RESULTS, "s28_A_chain_rows.jsonl"))
    out = {}
    for pdb in PDBS:
        t0 = time.time()
        r = A.run_chain_target(pdb, ["prod"])
        new, old = r["arms"]["prod"]["rmsd_chain"], rows[pdb]["arms"]["prod"]["rmsd_chain"]
        out[pdb] = dict(stored=old, reprojected=new, diff=new - old, bit_identical=(new == old), secs=time.time() - t0)
        print("  %s stored %.10f re-projected %.10f diff %+.3e %s (%.1fs)" % (pdb, old, new, new - old, "IDENTICAL" if new == old else "DIFFERS", time.time() - t0), flush=True)
    res = dict(targets=out, n=len(out), n_identical=int(sum(v["bit_identical"] for v in out.values())),
               max_abs_diff=max(abs(v["diff"]) for v in out.values()))
    path = os.path.join(A.RESULTS, "s28_A2_prodcheck.json")
    ST.save_atomic(path, res, module_file=__file__)
    print("identical on %d/%d, max |diff| %.3e; written %s" % (res["n_identical"], res["n"], res["max_abs_diff"], path))


if __name__ == "__main__":
    main()
