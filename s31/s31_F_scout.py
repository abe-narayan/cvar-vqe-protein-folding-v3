#!/usr/bin/env python
"""s31/s31_F_scout.py -- lane F scout: verify the reproduction gate, the operator identities,
and time the projection.  NO ENDPOINT NUMBER IS PRODUCED HERE.  Prereg: s31/PREREG_S31_F.md.
"""
from __future__ import annotations
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I  # noqa: E402

POOL_K, M, MIN_SEP = 500, 75, 2

for pdb in ("1ID6", "9KAR", "1A13"):
    t = [x for x in I.targets() if x["pdb"] == pdb][0]
    u = I.load_univ(pdb); rec = I.shipped_record(pdb)
    dg = I.distogram(pdb)
    pool = np.asarray(u["order"], int)[:POOL_K]
    W = np.asarray(u["W"], float)[pool]
    rr = np.asarray(u["rr"], float)[pool]
    nat = np.asarray(u["nat_ca"], float)
    ii, jj = I.pair_index(int(u["n"]), MIN_SEP)
    D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    sc = I.shipped_score(dg, D)
    order = np.argsort(sc, kind="stable")
    top = order[:M]
    sub = np.asarray(rec["sub"], int)
    gate = set(top.tolist()) == set(sub.tolist())
    Pt = I.pairwise_rmsd(W[top])
    b = I.medoid(Pt)
    Sup = I.superpose_batch(W[top], W[top[b]])
    C = Sup.mean(0)
    cl_avg = I.ca_rmsd(C, nat)
    cl_med = I.ca_rmsd(W[top[b]], nat)
    disp = float((Pt.sum() - np.trace(Pt)) / (M * (M - 1)))
    S = float(np.sqrt(((Sup - C) ** 2).sum(2).mean()))          # per-atom RMS spread
    s_b = float(I.ca_rmsd(Sup[b], C))
    t0 = time.time(); pa = I.project(C, t["seq"], t["fold"]); ta = time.time() - t0
    t0 = time.time(); pm = I.project(W[top[b]], t["seq"], t["fold"]); tm = time.time() - t0
    ch_avg = I.ca_rmsd(np.asarray(pa["ca"], float), nat)
    ch_med = I.ca_rmsd(np.asarray(pm["ca"], float), nat)
    # CA-CA bond geometry of the raw medoid (is it already a valid chain?)
    d_med = np.linalg.norm(np.diff(W[top[b]], axis=0), axis=1)
    d_avg = np.linalg.norm(np.diff(C, axis=0), axis=1)
    print("%s n=%2d gate=%s  rec.rmsd_avg %.4f  my cl_avg %.4f | rec.rmsd_fit %.4f my ch_avg %.4f"
          % (pdb, u["n"], gate, rec["rmsd_avg"], cl_avg, rec["rmsd_fit"], ch_avg))
    print("    cloud: AVG %.4f  MED %.4f  (med-avg %+.4f)   chain: AVG %.4f  MED %.4f (med-avg %+.4f)"
          % (cl_avg, cl_med, cl_med - cl_avg, ch_avg, ch_med, ch_med - ch_avg))
    print("    P(AVG) %+.4f  P(MED) %+.4f   DISP %.3f  S %.3f  s_b %.3f  B %.3f  sqrt(B^2+s_b^2) %.3f"
          % (ch_avg - cl_avg, ch_med - cl_med, disp, S, s_b, cl_avg, np.hypot(cl_avg, s_b)))
    print("    CA-CA bond: med mean %.3f sd %.4f | avg mean %.3f sd %.4f | proj %.2fs/%.2fs  n_distinct=%s"
          % (d_med.mean(), d_med.std(), d_avg.mean(), d_avg.std(), ta, tm, rec.get("n_distinct")))
