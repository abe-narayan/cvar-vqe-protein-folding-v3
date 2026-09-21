#!/usr/bin/env python
"""s31/s31_C_cache.py -- lane C's shared pool table.

One npz per target holding everything C1/C2/C3 need, so no downstream script re-loads a pool:

    W       (k, n, 3) float32   the shipped K=500 BLOSUM pool, universe frame
    order   (k,) int32          the DIS score order (production tie key), order[:75] == the top-75
    dis     (k,) float32        the shipped DIS score
    rr      (k,) float32        ORACLE Ca-RMSD of each member to the native.  ORACLE.
    nat     (n, 3) float32      the native CA trace.  ORACLE.
    top75   (75,) int32         the production record's retained set (sorted)
    n, k, fold

The ORACLE quantities live in the same file as the native-free ones deliberately: every consumer
labels its own arms, and hiding `rr` in a second file has never once stopped a leak in this project.

    python s31/s31_C_cache.py build [--limit N]
    python s31/s31_C_cache.py check
"""
from __future__ import annotations

import argparse
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

from s12 import instrument as I                    # noqa: E402

CACHE = os.path.join(HERE, "results", "s31_C_cache")


def path_for(pdb):
    return os.path.join(CACHE, "%s.npz" % pdb)


def build_one(pdb):
    from s29 import s29_O_ladder as O
    cand, dis, top, order, dg = O.load_pool(pdb)
    W = np.asarray(cand.W, np.float32)
    rr = np.asarray(cand.oracle_rr, np.float32)                 # ORACLE
    nat = np.asarray(cand.nat_ca, np.float32)                   # ORACLE
    assert set(np.asarray(order[:75], int).tolist()) == set(np.asarray(top, int).tolist()), \
        "%s: DIS order head != production top-75" % pdb
    os.makedirs(CACHE, exist_ok=True)
    tmp = path_for(pdb) + ".%d.tmp" % os.getpid()               # unique temp, never a shared name
    with open(tmp, "wb") as fh:                                 # handle, so numpy cannot rename it
        np.savez_compressed(fh, W=W, order=np.asarray(order, np.int32),
                            dis=np.asarray(dis, np.float32), rr=rr, nat=nat,
                            top75=np.asarray(top, np.int32),
                            n=np.int32(cand.n), k=np.int32(cand.k), fold=np.int32(cand.fold))
    os.replace(tmp, path_for(pdb))
    return int(cand.k), int(cand.n)


def load(pdb):
    z = np.load(path_for(pdb), allow_pickle=False)
    return {k: z[k] for k in z.files}


def all_pdbs():
    return [t["pdb"] for t in I.targets()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["build", "check"])
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    pdbs = all_pdbs()
    if a.limit:
        pdbs = pdbs[:a.limit]
    if a.phase == "check":
        miss = [p for p in pdbs if not os.path.exists(path_for(p))]
        print("cache: %d/%d present, %d missing" % (len(pdbs) - len(miss), len(pdbs), len(miss)))
        if miss:
            print("missing:", miss[:10])
        return
    t0 = time.time()
    for i, p in enumerate(pdbs):
        if os.path.exists(path_for(p)):
            continue
        k, n = build_one(p)
        if i % 10 == 0 or i == len(pdbs) - 1:
            print("[%3d/%3d] %s k=%d n=%2d  %.1fs" % (i + 1, len(pdbs), p, k, n, time.time() - t0),
                  flush=True)
    print("done %d targets in %.1fs" % (len(pdbs), time.time() - t0))


if __name__ == "__main__":
    main()
