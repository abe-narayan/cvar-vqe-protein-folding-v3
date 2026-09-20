#!/usr/bin/env python
"""s30/s30_P_prior.py -- S30 lane P, second measurement.

THE PRIOR ROUTE, PRICED ON THE PINNED BENCHMARK.  Lane L's S30-L18 relays S19's result that the
distogram's harmful mode is a per-target SEPARATION PROFILE worth -0.525 A ORACLE, with the best
native-free estimator recovering 19-30% of it.  **S19 predates the benchmark pinning**, which once
silently moved 13 targets and invalidated every fold model, so the coordinator asked for the number
to be re-checked on the pinned 126 before anything is built on it.

Every arm corrects the distogram by translating each pair's POSTERIOR by Delta_p, so the shipped
Bayes-risk score is evaluated unchanged at  risk_p(d - Delta_p).  The pipeline downstream is
production's: rescore the same 500-member BLOSUM pool, keep the top 75, coordinate-average in the
medoid frame.  The endpoint reported is the POINT-CLOUD Ca RMSD (production 3.0483 A) and is an
INTERMEDIATE; the built chain is 3.2041 A and the headline arms are projected to it.

ORACLE arms read the native distances and are labelled ORACLE in every sentence.  The deployable
arms read no native at inference and are fitted leave-fold-out on the pinned folds.

    python s30/s30_P_prior.py run
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
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s30_P_prior.json")
CACHE = os.path.join(RESULTS, "s30_P_prior_rows.npz")
POOL_K, M, MIN_SEP = 500, 75, 2
SEP_EDGES = [(2, 2), (3, 3), (4, 4), (5, 6), (7, 99)]      # the five-number profile, fixed bins


def bin_of(sep):
    b = np.zeros(len(sep), int)
    for k, (lo, hi) in enumerate(SEP_EDGES):
        b[(sep >= lo) & (sep <= hi)] = k
    return b


def score_shift(dg, D, delta):
    """The shipped Bayes-risk score with each pair's posterior translated by delta_p."""
    grid = dg["grid"]; risk = dg["risk"]
    g = np.clip(((np.asarray(D, float) - delta[None] - grid[0]) / 0.05).astype(np.int32),
                0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


def emit(W, sc, nat):
    """Production's readout: top-75 by score, coordinate average in the medoid frame."""
    top = np.argsort(sc, kind="stable")[:M]
    C, _ = I.coordinate_average(W[top])
    return float(I.ca_rmsd(C, nat)), C, top


def rows():
    out = []
    t0 = time.time()
    for k, t in enumerate(I.targets()):
        pdb = t["pdb"]; n = int(t["n"])
        u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
        pool = np.asarray(u["order"], int)[:POOL_K]
        W = np.asarray(u["W"], float)[pool]
        nat = np.asarray(u["nat_ca"], float)
        ii, jj = I.pair_index(n, MIN_SEP)
        exp = np.asarray(dg["expected"], float)
        d_nat = np.linalg.norm(nat[ii] - nat[jj], axis=1)                 # ORACLE
        D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        b = bin_of(jj - ii)
        err = exp - d_nat                                                  # ORACLE
        prof = np.array([err[b == q].mean() if (b == q).any() else 0.0 for q in range(5)])
        # per-residue additive (n numbers): least squares  err_p ~ r_i + r_j
        A = np.zeros((len(err), n))
        A[np.arange(len(err)), ii] += 1.0; A[np.arange(len(err)), jj] += 1.0
        r_res = np.linalg.lstsq(A, err, rcond=None)[0]
        # multiplicative stretch (1 number): exp*(1+a) ~ d_nat
        a_str = float((exp @ (d_nat - exp)) / max(exp @ exp, 1e-12))
        # native-free pool regression of the profile (what the pipeline already does)
        sub = np.asarray(rec["sub"], int)
        pool_err = (exp[None] - D[sub]).mean(0)
        prof_pool = np.array([pool_err[b == q].mean() if (b == q).any() else 0.0 for q in range(5)])
        out.append(dict(pdb=pdb, n=n, fold=int(t["fold"]), fail18=pdb in I.FAIL18,
                        prof=prof, prof_pool=prof_pool, a_str=a_str, off=float(err.mean()),
                        r_res=r_res, b=b, err=err, exp=exp, d_nat=d_nat, ii=ii, jj=jj,
                        pool_mean=float(np.asarray(u["rr"], float)[pool].mean()),
                        pool_best=float(np.asarray(u["rr"], float)[pool].min())))
        out[-1]["_W"] = W; out[-1]["_nat"] = nat; out[-1]["_D"] = D; out[-1]["_dg"] = dg
        if (k + 1) % 20 == 0:
            print("  rows %d/126 %.1fs" % (k + 1, time.time() - t0), flush=True)
    return out


def main():
    R = rows()
    folds = np.array([r["fold"] for r in R]); pdbs = [r["pdb"] for r in R]
    F = sorted(set(folds.tolist()))
    # leave-fold-out GLOBAL profile: one 5-vector shared by every target in the held-out fold
    gprof = {}
    for f in F:
        P = np.array([r["prof"] for r in R if r["fold"] != f])
        gprof[f] = P.mean(0)
    arms = {}

    def run(name, deltafn):
        vals, structs = [], []
        for r in R:
            d = deltafn(r)
            sc = score_shift(r["_dg"], r["_D"], -d)   # corrected prediction = exp - d
            v, C, top = emit(r["_W"], sc, r["_nat"])
            vals.append(v); structs.append(C)
        arms[name] = np.array(vals)
        return np.array(vals), structs

    base, base_s = run("PROD", lambda r: np.zeros_like(r["exp"]))
    run("ORACLE_FULL", lambda r: r["err"])
    run("ORACLE_SEPPROF5", lambda r: r["prof"][r["b"]])
    run("ORACLE_OFFSET1", lambda r: np.full_like(r["exp"], r["off"]))
    run("ORACLE_STRETCH1", lambda r: -r["a_str"] * r["exp"])
    run("ORACLE_PERRES", lambda r: r["r_res"][r["ii"]] + r["r_res"][r["jj"]])
    run("LFO_GLOBALPROF5", lambda r: gprof[r["fold"]][r["b"]])
    run("NF_POOLPROF5", lambda r: r["prof_pool"][r["b"]])

    nt = len(R)
    pm = np.array([r["pool_mean"] for r in R]); pb = np.array([r["pool_best"] for r in R])
    f18 = np.array([r["fail18"] for r in R])
    tA = np.zeros(nt, bool); tA[np.argsort(-pm)[:18]] = True
    tB = np.zeros(nt, bool); tB[np.argsort(-pb)[:18]] = True
    strata = {"all126": np.ones(nt, bool), "FAIL18_CIRCULAR": f18, "other108": ~f18,
              "tailA_poolmean": tA, "tailB_oraclebest": tB}

    res = {"meta": {"n": nt, "prod_cloud_mean": float(base.mean()),
                    "note": "POINT CLOUD Ca RMSD, an INTERMEDIATE; production built chain 3.2041 A",
                    "sep_bins": SEP_EDGES,
                    "lfo_global_profiles": {str(f): gprof[f].tolist() for f in F},
                    "oracle_profile_mean": np.array([r["prof"] for r in R]).mean(0).tolist(),
                    "oracle_profile_sd": np.array([r["prof"] for r in R]).std(0).tolist()},
           "arms": {}, "strata": {}}
    for k, v in arms.items():
        if k == "PROD":
            continue
        c = ST.compare(v, base, folds=folds, names=pdbs, label=k)
        res["arms"][k] = {"mean": float(v.mean()), "delta": c["effect"], "se": c["se"],
                          "mde": c["mde"], "x_mde": c["effect_over_mde"],
                          "ci95_fold": c["ci95_fold"], "ci95_iid": c["ci95_iid"],
                          "folds_same_sign": c["folds_same_sign"], "W": c["n_better"],
                          "L": c["n_worse"], "verdict": c["verdict"], "per_fold": c["per_fold"]}
        for sn, mk in strata.items():
            if sn == "all126":
                continue
            res["strata"].setdefault(k, {})[sn] = {
                "n": int(mk.sum()), "mean": float(v[mk].mean()),
                "delta": float((v - base)[mk].mean())}
    full = arms["ORACLE_FULL"] - base
    for k, v in arms.items():
        if k.startswith("ORACLE") and k != "ORACLE_FULL":
            res["arms"][k]["share_of_ORACLE_FULL_gap"] = float((v - base).mean() / full.mean())
    for k in ("LFO_GLOBALPROF5", "NF_POOLPROF5"):
        sp = arms["ORACLE_SEPPROF5"] - base
        res["arms"][k]["share_of_ORACLE_SEPPROF5"] = float((arms[k] - base).mean() / sp.mean())
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for k, a in res["arms"].items():
        print("%-18s %.4f  delta %+.4f  %.2fxMDE  fold95 [%+.4f,%+.4f]  %dW/%dL  ss%d  %s"
              % (k, a["mean"], a["delta"], abs(a["x_mde"]), a["ci95_fold"][0], a["ci95_fold"][1],
                 a["W"], a["L"], a["folds_same_sign"], a["verdict"]))
    print("production point cloud %.4f" % base.mean())
    print("wrote", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("stage", nargs="?", default="run")
    ap.parse_args(); main()
