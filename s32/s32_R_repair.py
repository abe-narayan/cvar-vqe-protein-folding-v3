#!/usr/bin/env python
"""S32 LANE R -- R1 (is the +0.1622 projection price recoverable?) and R4 (can the readout
land in the cheap projection regime?), every arm projected IN THE SAME PROCESS.

PREREG: s32/PREREG_S32_R.md at 02754f5a.

ARMS (all use `core.project.lam_path`, i.e. PRODUCTION'S OWN CODE PATH, ramah/lam=0.3/
multi=True/maxiter=300/grad=exact; the only thing that varies is the cloud handed to it and
whether one extra start is offered):

    PROD            lam_path(C,               extra=None)   <- production, the baseline
    MEDOID_EXTRA    lam_path(C,               extra=medoid torsions)      DEPLOYABLE
    SCALE_NF        lam_path(s_nf * C,        extra=None)                 DEPLOYABLE
    SCALE_NF_MED    lam_path(s_nf * C,        extra=medoid torsions)      DEPLOYABLE
    SCALE_GRID_*    lam_path(s * C, extra=None) for s on a fixed grid
                    -> SCALE_ORACLE = per-target argmin on NATIVE RMSD.  ORACLE / NOT
                       DEPLOYABLE, prices the ceiling of the scale axis only.

`s_nf = ideal_virtual_bond / mean_virtual_bond(C)` is computed from the cloud's own geometry.
No native, no label, no fold model.  The scaling is about the ORIGIN-CENTRED cloud so it is
a pure dilation of shape.

CONTRACT RULE 15: every arm carries its emitted chain's virtual-bond mean/sd and its
displacement from the production chain and from the cloud it was projected from, from row one.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                          # noqa: E402
from core import project as pj                                           # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
S29_STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")
os.makedirs(RESULTS, exist_ok=True)

LAM, MAXITER, GRAD = 0.3, 300, "exact"
SCALE_GRID = (0.94, 0.97, 1.00, 1.03, 1.06, 1.09, 1.12, 1.15)


def ideal_bond():
    """The builder's virtual CA-CA distance.  A constant of the ideal geometry, measured
    from the shipped builder rather than quoted."""
    ca = np.asarray(pj.build_ca_exact(np.full((1, 8), -1.0), np.full((1, 8), 2.0)), float)[0]
    return float(np.linalg.norm(ca[1:] - ca[:-1], axis=-1).mean())


IDEAL = ideal_bond()


def dilate(C, s):
    m = C.mean(0)
    return m + s * (C - m)


def top75(u, pdb):
    pool = I.pool_idx(u)
    sub = np.asarray(I.shipped_record(pdb)["sub"], int)
    idx = pool[sub]
    W = np.asarray(u["W"], float)[idx].astype(np.float32).astype(float)
    P = I.pairwise_rmsd(W).astype(np.float32).astype(float)
    b = int(np.argmin(P.sum(1) / max(P.shape[0] - 1, 1)))
    return idx, b


def arm(C, seq, fold, extra, nat, Cproj, prod_ca):
    pen = pj.make_penalty("ramah", seq, int(fold))
    path = pj.lam_path(np.asarray(C, float), pen, (0.0, LAM), extra=extra,
                       maxiter=MAXITER, multi=True, grad=GRAD)
    ca = np.asarray(path[LAM][0], float)
    vb = np.linalg.norm(ca[1:] - ca[:-1], axis=-1)
    out = dict(rmsd=float(I.ca_rmsd(ca, nat)),
               fit_rmsd=float(I.ca_rmsd(np.asarray(path[0.0][0], float), nat)),
               obj=float(path[LAM][3]), d_to_C=float(path[LAM][5]),
               vbond_mean=float(vb.mean()), vbond_sd=float(vb.std()),
               d_to_cloud=float(I.ca_rmsd(ca, Cproj)))
    out["d_to_prod"] = float(I.ca_rmsd(ca, prod_ca)) if prod_ca is not None else 0.0
    return out, ca


def one(pdb, verbose=True):
    t0 = time.time()
    with np.load(os.path.join(S29_STRUCTS, f"{pdb}.npz"), allow_pickle=True) as z:
        C = np.asarray(z["prod"], float)
        seq, fold, n = str(z["seq"]), int(z["fold"]), int(z["n"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)                   # ORACLE: scoring only
    idx, mpos = top75(u, pdb)
    ph_m = np.asarray(u["PHI"], float)[idx[mpos]].copy()
    ps_m = np.asarray(u["PSI"], float)[idx[mpos]].copy()

    vb_cloud = float(np.linalg.norm(C[1:] - C[:-1], axis=-1).mean())
    s_nf = IDEAL / vb_cloud
    Cs = dilate(C, s_nf)

    out = dict(pdb=pdb, n=n, fold=fold, ideal_bond=IDEAL,
               cloud_vbond_mean=vb_cloud,
               cloud_vbond_sd=float(np.linalg.norm(C[1:] - C[:-1], axis=-1).std()),
               s_nf=s_nf, medoid_pos=mpos,
               cloud_rmsd=float(I.ca_rmsd(C, nat)),
               cloud_scaled_rmsd=float(I.ca_rmsd(Cs, nat)))

    a_prod, prod_ca = arm(C, seq, fold, None, nat, C, None)
    out["PROD"] = a_prod
    out["MEDOID_EXTRA"] = arm(C, seq, fold, (ph_m, ps_m), nat, C, prod_ca)[0]
    out["SCALE_NF"] = arm(Cs, seq, fold, None, nat, Cs, prod_ca)[0]
    out["SCALE_NF_MED"] = arm(Cs, seq, fold, (ph_m, ps_m), nat, Cs, prod_ca)[0]
    grid = {}
    for s in SCALE_GRID:
        Cg = dilate(C, s)
        grid["%.2f" % s] = dict(arm(Cg, seq, fold, None, nat, Cg, prod_ca)[0],
                                cloud_rmsd=float(I.ca_rmsd(Cg, nat)))
    out["SCALE_GRID"] = grid
    out["secs"] = round(time.time() - t0, 2)
    if verbose:
        print("%-6s n=%2d s_nf=%.4f PROD %.4f MED %.4f SCALE %.4f SCALE+MED %.4f "
              "gridmin %.4f  %.0fs"
              % (pdb, n, s_nf, out["PROD"]["rmsd"], out["MEDOID_EXTRA"]["rmsd"],
                 out["SCALE_NF"]["rmsd"], out["SCALE_NF_MED"]["rmsd"],
                 min(v["rmsd"] for v in grid.values()), out["secs"]), flush=True)
    return out


def rows_path(shard):
    return os.path.join(RESULTS, f"s32_R_repair_shard{int(shard)}.jsonl")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    pdbs = [t["pdb"] for t in I.targets()]
    mine = [p for k, p in enumerate(pdbs) if k % a.nshard == a.shard]
    if a.limit:
        mine = mine[:a.limit]
    path = rows_path(a.shard)
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["pdb"])
                except Exception:                                        # noqa: BLE001
                    pass
    for p in mine:
        if p in done:
            continue
        r = one(p)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
    print("shard %d done (%d targets) -> %s" % (a.shard, len(mine), path), flush=True)


if __name__ == "__main__":
    main()
