#!/usr/bin/env python
"""s30/s30_Q_rankcheck.py -- lane Q's INDEPENDENT verification of lane T's stable-rank number.

Deliberately a different data path from `s30/s30_T_bits.py`: T reads `s8/generate_univ/<pdb>.npz`
and poses on its own medoid; this reads the DEPLOYED pool through `s29_O_ladder.load_pool` and the
deployed `s28_A_amp.Frame`, which is the frame the readout actually uses. If the two agree, the
number is robust to the path; if they disagree, the disagreement is the finding.

Also reports the stable rank in COORDINATE space, because the sparse/weighted readout and any
second-moment (quadric) construction act on coordinates, not on the pair-distance map.

    python s30/s30_Q_rankcheck.py [--limit N]
"""
from __future__ import annotations

import argparse
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
from s29 import s29_O_ladder as O          # noqa: E402

OUT = os.path.join(HERE, "results", "s30_Q_rankcheck.json")


def stable_rank(A):
    """||A||_F^2 / ||A||_2^2 on the column-centred A, plus the top-eigenvalue share and k90."""
    Ac = np.asarray(A, float)
    Ac = Ac - Ac.mean(0, keepdims=True)
    s = np.linalg.svd(Ac, compute_uv=False)
    lam = s ** 2
    if lam[0] <= 0:
        return 1.0, 1.0, 1
    tot = lam.sum()
    k90 = int(np.searchsorted(np.cumsum(lam) / tot, 0.90) + 1)
    return float(tot / lam[0]), float(lam[0] / tot), k90


def main(limit=None):
    pdbs = O.all_pdbs()
    if limit:
        pdbs = pdbs[:int(limit)]
    rows = []
    for pdb in pdbs:
        cand, dis, top, order, dg = O.load_pool(pdb)
        W = cand.W
        n = cand.n
        from s27 import s28_A_amp as A
        frame = A.Frame(W, np.sort(top))          # the DEPLOYED frame
        # coordinate feature matrix, in the frame the readout uses
        rs_c, top_c, k90_c = stable_rank(frame.Wf)
        # pair-distance feature matrix, min_sep = 2 (lane T's convention)
        i, j = np.triu_indices(n, k=2)
        D = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=2)
        rs_d, top_d, k90_d = stable_rank(D)
        rows.append(dict(pdb=pdb, n=n, k=int(cand.k), npairs=int(len(i)),
                         rs_coord=rs_c, top_coord=top_c, k90_coord=k90_c,
                         rs_dist=rs_d, top_dist=top_d, k90_dist=k90_d))
        print("  %-6s n=%2d  r_stable(dist) %.3f (PC1 %.1f%%, k90 %d) | "
              "r_stable(coord) %.3f (PC1 %.1f%%, k90 %d)"
              % (pdb, n, rs_d, 100 * top_d, k90_d, rs_c, 100 * top_c, k90_c), flush=True)

    rd = np.array([r["rs_dist"] for r in rows])
    rc = np.array([r["rs_coord"] for r in rows])
    out = dict(n=len(rows), rows=rows,
               rs_dist_mean=float(rd.mean()), rs_dist_median=float(np.median(rd)),
               rs_dist_max=float(rd.max()), rs_dist_min=float(rd.min()),
               rs_coord_mean=float(rc.mean()), rs_coord_median=float(np.median(rc)),
               rs_coord_max=float(rc.max()), rs_coord_min=float(rc.min()),
               pc1_dist_mean=float(np.mean([r["top_dist"] for r in rows])),
               k90_dist_mean=float(np.mean([r["k90_dist"] for r in rows])),
               pc1_coord_mean=float(np.mean([r["top_coord"] for r in rows])),
               k90_coord_mean=float(np.mean([r["k90_coord"] for r in rows])),
               provenance=ST.provenance(__file__))
    print("\n== lane Q independent rank check, n = %d pools, DEPLOYED pool + DEPLOYED frame" % len(rows))
    print("   PAIR-DISTANCE  r_stable mean %.3f  median %.3f  max %.3f   PC1 %.1f%%  k90 %.1f"
          % (out["rs_dist_mean"], out["rs_dist_median"], out["rs_dist_max"],
             100 * out["pc1_dist_mean"], out["k90_dist_mean"]))
    print("   COORDINATE     r_stable mean %.3f  median %.3f  max %.3f   PC1 %.1f%%  k90 %.1f"
          % (out["rs_coord_mean"], out["rs_coord_median"], out["rs_coord_max"],
             100 * out["pc1_coord_mean"], out["k90_coord_mean"]))
    print("   lane T reported 1.859 / 1.865 / 2.72 (dist) and 3.404 (coord) from s8/generate_univ.")
    ST.save_atomic(OUT, out, module_file=__file__)
    print("wrote", OUT)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    main(a.limit)
