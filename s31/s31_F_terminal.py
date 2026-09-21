#!/usr/bin/env python
"""s31/s31_F_terminal.py -- S31 LANE F: the terminal operator, three arms on the BUILT CHAIN.

Pre-registered in `s31/PREREG_S31_F.md` (commits 8a14edea + 49ee7c92, both before this file).

Three terminal operators over the SAME shipped top-75 of the SAME 500-member BLOSUM pool:

    AVG     production: superpose on the unweighted consensus medoid, uniform coordinate mean
    MED     the quantum stage's readout: return the consensus medoid MEMBER (a real backbone)
    AVG_RG  AVG rescaled about its own centroid to the 75 members' mean Rg  (native-free)

all passed through the SAME stage-3b projection `I.project` in the SAME process, so contract
rule 31 (unpinned projection seed, 0.0107 A instrument spread) cannot decide a comparison.

REPRODUCTION GATE: the recomputed shipped top-75 must equal the production record's `sub`
set-wise, and the recomputed AVG cloud RMSD must match `rec["rmsd_avg"]` to 1e-3, on every
target, or the run aborts.

ORACLE: every quantity derived from `rr` or `nat_ca` is an ORACLE diagnostic.  No deployable
parameter in this file is chosen on the native.

    python s31/s31_F_terminal.py run        # per-target, checkpointed to the jsonl
    python s31/s31_F_terminal.py analyse
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

from s12 import instrument as I            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s31_F_terminal_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_F_terminal.json")

POOL_K, M, MIN_SEP, SEED = 500, 75, 2, 31006


def rg(X):
    X = np.asarray(X, float)
    return float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))


def bond(X):
    d = np.linalg.norm(np.diff(np.asarray(X, float), axis=0), axis=1)
    return float(d.mean()), float(d.std())


def sepprof(X, n):
    """Mean CA-CA distance at each sequence separation s = 1..n-1.  The corrected averaging
    mechanism (s15/coord_FINDINGS.md:914-921) is a SEPARATION-DEPENDENT shape distortion --
    short contracted, long expanded, crossing unity near |i-j| = 8 -- not a uniform scale."""
    X = np.asarray(X, float)
    d = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    i, j = np.triu_indices(n, 1)
    s = j - i
    return [float(d[i[s == q], j[s == q]].mean()) for q in range(1, n)]


def sepprof_batch(W, n):
    """The same profile for a STACK of structures, averaged over members (native-free)."""
    W = np.asarray(W, float)
    d = np.linalg.norm(W[:, :, None, :] - W[:, None, :, :], axis=3)
    i, j = np.triu_indices(n, 1)
    s = j - i
    return [float(d[:, i[s == q], j[s == q]].mean()) for q in range(1, n)]


def mds(Dm):
    """Classical MDS: the best 3-D point set for a target distance matrix (eigenvalues clipped)."""
    Dm = np.asarray(Dm, float)
    n = Dm.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (Dm ** 2) @ J
    w, V = np.linalg.eigh(B)
    o = np.argsort(-w)[:3]
    return V[:, o] * np.sqrt(np.maximum(w[o], 0.0))[None, :]


def one_target(t):
    pdb = t["pdb"]
    u = I.load_univ(pdb)
    rec = I.shipped_record(pdb)
    dg = I.distogram(pdb)
    n = int(u["n"])
    pool = np.asarray(u["order"], int)[:POOL_K]
    W = np.asarray(u["W"], float)[pool]
    rr = np.asarray(u["rr"], float)[pool]                    # ORACLE: per-member RMSD to native
    nat = np.asarray(u["nat_ca"], float)
    ii, jj = I.pair_index(n, MIN_SEP)
    D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    sc = I.shipped_score(dg, D)
    order = np.argsort(sc, kind="stable")
    top = order[:M]
    if set(top.tolist()) != set(np.asarray(rec["sub"], int).tolist()):
        raise RuntimeError("REPRODUCTION GATE FAILED (top-75 != production sub) on %s" % pdb)

    Pt = I.pairwise_rmsd(W[top])
    b = int(I.medoid(Pt))
    Sup = I.superpose_batch(W[top], W[top[b]])
    C = Sup.mean(0)
    cl_avg = I.ca_rmsd(C, nat)
    if abs(cl_avg - float(rec["rmsd_avg"])) > 1e-3:
        raise RuntimeError("REPRODUCTION GATE FAILED (cloud %.6f vs rec %.6f) on %s"
                           % (cl_avg, rec["rmsd_avg"], pdb))
    Xmed = W[top[b]]
    cl_med = I.ca_rmsd(Xmed, nat)

    # ---- AVG_RG: rescale the average about its own centroid to the members' mean Rg (native-free)
    rg_members = np.sqrt(((W[top] - W[top].mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    rg_target = float(rg_members.mean())
    rg_C = rg(C)
    lam = rg_target / max(rg_C, 1e-9)
    Crg = (C - C.mean(0)) * lam + C.mean(0)
    cl_rg = I.ca_rmsd(Crg, nat)

    # ---- AVG_SEP: rescale the average's pair distances PER SEQUENCE SEPARATION to the members'
    # own per-separation means, then re-embed.  NATIVE-FREE.  This is the correction matched to
    # the CORRECTED mechanism (separation-dependent shape distortion), where AVG_RG is the
    # uniform-scale correction matched to the WITHDRAWN 25.8% contraction story.
    p_mem = np.array(sepprof_batch(W[top], n), float)        # native-free target profile
    p_avg = np.array(sepprof(C, n), float)
    f = p_mem / np.maximum(p_avg, 1e-9)
    Dc = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=2)
    sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    F = np.ones((n, n))
    nz = sep > 0
    F[nz] = f[sep[nz] - 1]
    Csep = mds(Dc * F)
    cl_sep = I.ca_rmsd(Csep, nat)

    # ---- stage 3b, identical for all arms
    out = {}
    for nm, X in (("AVG", C), ("MED", Xmed), ("AVG_RG", Crg), ("AVG_SEP", Csep)):
        X = np.asarray(X, float)
        pr = I.project(X, t["seq"], t["fold"])
        ca = np.asarray(pr["ca"], float)
        # MOVE = how far stage 3b has to drag the input to make it a chain.  NATIVE-FREE:
        # both arguments are available at inference.  See the prereg's second addendum.
        out[nm] = dict(chain=float(I.ca_rmsd(ca, nat)), rg_chain=rg(ca),
                       bond_chain=bond(ca), move=float(I.ca_rmsd(ca, X)),
                       prof=sepprof(ca, n))

    disp = float((Pt.sum() - np.trace(Pt)) / (M * (M - 1)))
    disp_med = float(np.median(Pt[np.triu_indices(M, 1)]))
    S = float(np.sqrt(((Sup - C) ** 2).sum(2).mean()))
    s_b = float(I.ca_rmsd(Sup[b], C))
    rank_best = int(np.where(order == int(np.argmin(rr)))[0][0])
    rank_best75 = int(np.where(order == int(top[np.argmin(rr[top])]))[0][0])

    row = dict(
        pdb=pdb, n=n, fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
        # --- ORACLE pool statistics
        pool_mean=float(rr.mean()), pool_best=float(rr.min()),
        set_mean=float(rr[top].mean()), set_best=float(rr[top].min()),
        rank_best_in_pool=rank_best, rank_best_in_top75=rank_best75,
        frac_top75_better_than_avg=float((rr[top] < cl_avg).mean()),
        n_top75_under3=int((rr[top] < 3.0).sum()),
        n_distinct=int(rec.get("n_distinct", -1)),
        # --- native-free set statistics
        DISP=disp, DISP_median=disp_med, S=S, s_b=s_b,
        rg_members_mean=rg_target, rg_members_sd=float(rg_members.std(ddof=1)),
        rg_scale_lambda=float(lam),
        # --- clouds
        cloud_AVG=cl_avg, cloud_MED=cl_med, cloud_AVG_RG=cl_rg, cloud_AVG_SEP=cl_sep,
        rg_nat=rg(nat), rg_cloud_AVG=rg_C, rg_cloud_MED=rg(Xmed), rg_cloud_AVG_RG=rg(Crg),
        bond_cloud_AVG=bond(C), bond_cloud_MED=bond(Xmed), bond_cloud_AVG_RG=bond(Crg),
        # --- separation profiles (mean CA-CA distance at each s = 1..n-1)
        prof_nat=sepprof(nat, n), prof_members=list(p_mem), prof_cloud_AVG=list(p_avg),
        prof_cloud_MED=sepprof(Xmed, n), prof_cloud_AVG_SEP=sepprof(Csep, n),
        prof_chain_AVG=out["AVG"]["prof"], prof_chain_MED=out["MED"]["prof"],
        prof_chain_AVG_SEP=out["AVG_SEP"]["prof"],
        # --- chains
        chain_AVG=out["AVG"]["chain"], chain_MED=out["MED"]["chain"],
        chain_AVG_RG=out["AVG_RG"]["chain"], chain_AVG_SEP=out["AVG_SEP"]["chain"],
        rg_chain_AVG=out["AVG"]["rg_chain"], rg_chain_MED=out["MED"]["rg_chain"],
        rg_chain_AVG_RG=out["AVG_RG"]["rg_chain"],
        move_AVG=out["AVG"]["move"], move_MED=out["MED"]["move"],
        move_AVG_RG=out["AVG_RG"]["move"], move_AVG_SEP=out["AVG_SEP"]["move"],
        bond_chain_AVG=out["AVG"]["bond_chain"], bond_chain_MED=out["MED"]["bond_chain"],
        rec_fit=float(rec["rmsd_fit"]), rec_avg=float(rec["rmsd_avg"]),
        prod_chain_record=float(rec["rmsd_fit"]),
    )
    row["P_AVG"] = row["chain_AVG"] - row["cloud_AVG"]
    row["P_MED"] = row["chain_MED"] - row["cloud_MED"]
    row["P_AVG_RG"] = row["chain_AVG_RG"] - row["cloud_AVG_RG"]
    row["P_AVG_SEP"] = row["chain_AVG_SEP"] - row["cloud_AVG_SEP"]
    return row


def run():
    done = set()
    if os.path.exists(ROWS):
        with open(ROWS) as fh:
            for ln in fh:
                if ln.strip():
                    done.add(json.loads(ln)["pdb"])
    ts = I.targets()
    t0 = time.time()
    for k, t in enumerate(ts):
        if t["pdb"] in done:
            continue
        r = one_target(t)
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        print("[%3d/%3d %5.0fs] %s n=%2d DISP %5.2f | cloud A %6.3f M %6.3f R %6.3f S %6.3f | "
              "chain A %6.3f M %6.3f R %6.3f S %6.3f | P_A %+6.3f P_M %+6.3f mv_A %5.3f"
              % (k + 1, len(ts), time.time() - t0, t["pdb"], r["n"], r["DISP"],
                 r["cloud_AVG"], r["cloud_MED"], r["cloud_AVG_RG"], r["cloud_AVG_SEP"],
                 r["chain_AVG"], r["chain_MED"], r["chain_AVG_RG"], r["chain_AVG_SEP"],
                 r["P_AVG"], r["P_MED"], r["move_AVG"]), flush=True)
    print("run complete:", len(ts), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["run"])
    a = ap.parse_args()
    run()
