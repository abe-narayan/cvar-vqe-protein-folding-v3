#!/usr/bin/env python
"""s30/s30_P_chain.py -- S30 lane P: the two headline ORACLE arms on the BUILT CHAIN.

Contract rule 1: built-chain Ca RMSD is the endpoint; point-cloud numbers are intermediates.
`s30_P_prior.py` reports the point cloud (production 3.0483 A).  This projects three arms through
production's stage-3b projection (`I.project`, verified to reproduce the record's `rmsd_fit` to
0.004 A) so the headline carries the endpoint.

    python s30/s30_P_chain.py
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s30.s30_P_prior import bin_of, score_shift, emit, POOL_K, M, MIN_SEP   # noqa: E402

OUT = os.path.join(HERE, "results", "s30_P_chain.json")
ROWS = os.path.join(HERE, "results", "s30_P_chain_rows.jsonl")


def main():
    F = [0, 1, 2, 3, 4]
    tg = I.targets()
    # first pass: the ORACLE separation profiles (needed for the leave-fold-out magnitude)
    prof, meta = {}, {}
    for t in tg:
        pdb = t["pdb"]; n = int(t["n"])
        u = I.load_univ(pdb); dg = I.distogram(pdb)
        ii, jj = I.pair_index(n, MIN_SEP)
        nat = np.asarray(u["nat_ca"], float)
        d_nat = np.linalg.norm(nat[ii] - nat[jj], axis=1)              # ORACLE
        err = np.asarray(dg["expected"], float) - d_nat                # ORACLE
        b = bin_of(jj - ii)
        prof[pdb] = np.array([err[b == q].mean() if (b == q).any() else 0.0 for q in range(5)])
        meta[pdb] = (b, ii, jj, int(t["fold"]))
    gmag = {f: np.abs(np.array([prof[t["pdb"]] for t in tg if t["fold"] != f])).mean(0) for f in F}

    done = set()
    if os.path.exists(ROWS):
        for ln in open(ROWS):
            done.add(json.loads(ln)["pdb"])
    fh = open(ROWS, "a")
    t0 = time.time()
    for k, t in enumerate(tg):
        pdb = t["pdb"]
        if pdb in done:
            continue
        u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
        b, ii, jj, fold = meta[pdb]
        pool = np.asarray(u["order"], int)[:POOL_K]
        W = np.asarray(u["W"], float)[pool]
        nat = np.asarray(u["nat_ca"], float)
        D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        row = dict(pdb=pdb, fold=fold, n=int(t["n"]), fail18=pdb in I.FAIL18,
                   pool_mean=float(np.asarray(u["rr"], float)[pool].mean()),
                   pool_best=float(np.asarray(u["rr"], float)[pool].min()),
                   rec_fit=float(rec["rmsd_fit"]))
        deltas = {"PROD": np.zeros(len(ii)),
                  "ORACLE_SEPPROF5": prof[pdb][b],
                  "ORACLEsign_LFOmag": np.sign(prof[pdb])[b] * gmag[fold][b]}
        for nm, d in deltas.items():
            sc = score_shift(dg, D, -d)
            v, C, _ = emit(W, sc, nat)
            pr = I.project(C, t["seq"], fold)
            row[nm + "_cloud"] = v
            row[nm + "_chain"] = float(I.ca_rmsd(np.asarray(pr["ca"], float), nat))
        fh.write(json.dumps(row) + "\n"); fh.flush()
        if (k + 1) % 10 == 0:
            print("  %d/%d  %.1fs" % (k + 1, len(tg), time.time() - t0), flush=True)
    fh.close()

    R = [json.loads(l) for l in open(ROWS)]
    R = {r["pdb"]: r for r in R}
    R = [R[t["pdb"]] for t in tg if t["pdb"] in R]
    folds = np.array([r["fold"] for r in R]); pdbs = [r["pdb"] for r in R]
    nt = len(R)
    pm = np.array([r["pool_mean"] for r in R]); pb = np.array([r["pool_best"] for r in R])
    tA = np.zeros(nt, bool); tA[np.argsort(-pm)[:18]] = True
    tB = np.zeros(nt, bool); tB[np.argsort(-pb)[:18]] = True
    f18 = np.array([r["fail18"] for r in R])
    strata = {"FAIL18_CIRCULAR": f18, "other108": ~f18, "tailA_poolmean": tA,
              "tailB_oraclebest": tB}
    base = np.array([r["PROD_chain"] for r in R])
    res = {"n": nt, "PROD_chain_mean": float(base.mean()),
           "PROD_cloud_mean": float(np.mean([r["PROD_cloud"] for r in R])),
           "record_fit_mean": float(np.mean([r["rec_fit"] for r in R])),
           "arms": {}, "strata": {}}
    for nm in ("ORACLE_SEPPROF5", "ORACLEsign_LFOmag"):
        v = np.array([r[nm + "_chain"] for r in R])
        c = ST.compare(v, base, folds=folds, names=pdbs, label=nm + "|chain")
        res["arms"][nm] = {"chain_mean": float(v.mean()), "delta": c["effect"], "se": c["se"],
                           "mde": c["mde"], "x_mde": c["effect_over_mde"],
                           "ci95_fold": c["ci95_fold"], "folds_same_sign": c["folds_same_sign"],
                           "W": c["n_better"], "L": c["n_worse"], "verdict": c["verdict"],
                           "cloud_delta": float(np.mean([r[nm + "_cloud"] - r["PROD_cloud"]
                                                         for r in R]))}
        for sn, mk in strata.items():
            res["strata"].setdefault(nm, {})[sn] = {
                "n": int(mk.sum()), "chain_mean": float(v[mk].mean()),
                "delta": float((v - base)[mk].mean())}
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(json.dumps(res["arms"], indent=1))
    print("PROD chain %.4f (record %.4f)  cloud %.4f"
          % (res["PROD_chain_mean"], res["record_fit_mean"], res["PROD_cloud_mean"]))


if __name__ == "__main__":
    main()
