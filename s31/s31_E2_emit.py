#!/usr/bin/env python
"""s31/s31_E2_emit.py -- S31 lane E: apply every correction through production's readout.

For each target and each arm in `s31_E2_deltas.npz`:
    score = score_shift(dg, D, -delta)   ->   top-75   ->   coordinate average in the medoid frame
    CA point-cloud RMSD (production 3.0483 A)        and
    I.project(...)  ->  BUILT-CHAIN CA RMSD (production 3.2105 A, S29 prod row)
Both bases are recorded for every arm; the chain is the endpoint and the cloud is an intermediate.

Shardable:  python s31/s31_E2_emit.py --shard 0 --nshards 5
Rows append to s31/results/s31_E2_rows.jsonl and the script resumes: a (pdb) already present is
skipped, so a killed shard can simply be relaunched.
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                                  # noqa: E402
from s30.s30_P_prior import score_shift, emit, POOL_K, MIN_SEP   # noqa: E402
from s31 import s31_E_lib as L                                   # noqa: E402

DELTAS = os.path.join(HERE, "results", "s31_E2_deltas.npz")
ROWS = os.path.join(HERE, "results", "s31_E2_rows.jsonl")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--rows", default=ROWS)
    a = ap.parse_args()

    z = np.load(DELTAS, allow_pickle=False)
    off = z["off"]; names = [str(s) for s in z["names"]]
    D_ = {k: z["delta_" + k] for k in names}
    tg = I.targets()
    assert len(tg) == len(off) - 1, "target count %d != cache %d" % (len(tg), len(off) - 1)

    rows_path = a.rows if a.nshards == 1 else a.rows.replace(".jsonl", ".s%d.jsonl" % a.shard)
    done = set()
    if os.path.exists(rows_path):
        for ln in open(rows_path):
            try:
                done.add(json.loads(ln)["pdb"])
            except Exception:
                pass
    fh = open(rows_path, "a")
    t0 = time.time(); k_done = 0
    for k, t in enumerate(tg):
        if k % a.nshards != a.shard:
            continue
        pdb = t["pdb"]
        if pdb in done:
            continue
        n = int(t["n"]); fold = int(t["fold"])
        u = I.load_univ(pdb); dg = I.distogram(pdb)
        ii, jj = I.pair_index(n, MIN_SEP)
        pool = np.asarray(u["order"], int)[:POOL_K]
        W = np.asarray(u["W"], float)[pool]
        nat = np.asarray(u["nat_ca"], float)
        Dm = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        s, e = int(off[k]), int(off[k + 1])
        assert e - s == len(ii), "pair count mismatch on %s" % pdb
        row = {"pdb": pdb, "fold": fold, "n": n, "fail18": pdb in I.FAIL18}
        for nm in names:
            sc = score_shift(dg, Dm, -D_[nm][s:e])
            v, C, _ = emit(W, sc, nat)
            pr = I.project(C, t["seq"], fold)
            row[nm + "_cloud"] = float(v)
            row[nm + "_chain"] = float(I.ca_rmsd(np.asarray(pr["ca"], float), nat))
        fh.write(json.dumps(row) + "\n"); fh.flush()
        k_done += 1
        print("  shard %d: %s (%d done) %.1fs" % (a.shard, pdb, k_done, time.time() - t0),
              flush=True)
    fh.close()
    print("shard %d finished, %d targets, %.1fs" % (a.shard, k_done, time.time() - t0))


if __name__ == "__main__":
    main()
