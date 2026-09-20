#!/usr/bin/env python
"""s30/s30_F_width.py -- S30 lane F, experiment F2: filter width vs averaging width.

Pre-registration `s30/PREREG_S30_F2.md` (committed e626a88e, before this file existed).

In production the top-75 is BOTH the filter and the averaging set.  This separates them:

    filter width  k : keep the top k of the 500-member BLOSUM pool by the shipped
                      distogram Bayes-risk score (s12.instrument.shipped_score)
    averaging width m : uniformly coordinate-average a RANDOM m-subset of those k

k = m = 75 reproduces production; k = 500, m = 75 is the "no filter" arm.  Holding m = 75 and
sweeping k moves the FILTER alone -- S29's flat m sweep was the AVERAGING width only.

REPRODUCTION GATE: the recomputed score's top-75 must equal the production record's `sub`
set-wise on all 126 targets, or the run aborts.  (Verified bit-exact on 1A13 and 9KAR first.)

ORACLE / DEPLOYABLE: `rr`-derived set statistics and every RMSD-to-native here are ORACLE
diagnostics; no deployable parameter is chosen on them in this file.  All RMSDs in phases 1-2
are POINT CLOUD and are labelled so; the built chain is phase 3 (`s30_F_width_chain.py`).

    python s30/s30_F_width.py cloud      # phases 1+2, per-target checkpointed
    python s30/s30_F_width.py analyse
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from zlib import crc32

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
ROWS = os.path.join(RESULTS, "s30_F_width_rows.jsonl")
OUT = os.path.join(RESULTS, "s30_F_width_cloud.json")

K_GRID = (75, 100, 128, 150, 200, 250, 300, 400, 500)      # prereg, fixed
M_GRID = (25, 50, 75, 100, 150)                            # prereg, fixed
N_SEED = 8                                                 # prereg
SEED = 30003                                               # prereg
POOL_K = 500
MIN_SEP = 2


def score_order(pdb, u, rec):
    """The shipped Bayes-risk score order over the 500-member pool, with the reproduction gate."""
    dg = I.distogram(pdb)
    pool = np.asarray(u["order"], int)[:POOL_K]
    W = np.asarray(u["W"], float)[pool]
    ii, jj = I.pair_index(int(u["n"]), min_sep=MIN_SEP)
    D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    sc = I.shipped_score(dg, D)
    order = np.argsort(sc, kind="stable")
    sub = np.asarray(rec["sub"], int)
    if set(order[:75].tolist()) != set(sub.tolist()):
        raise RuntimeError("REPRODUCTION GATE FAILED on %s: recomputed top-75 != production sub" % pdb)
    return order, W, np.asarray(u["rr"], float)[pool], sc


def one_target(pdb, t):
    u = I.load_univ(pdb)
    rec = I.shipped_record(pdb)
    order, W, rr_pool, sc = score_order(pdb, u, rec)
    nat = np.asarray(u["nat_ca"], float)
    P = I.pairwise_rmsd(W)                       # 500 x 500, computed once, sliced per cell
    rng = np.random.default_rng(SEED + crc32(pdb.encode()))   # deterministic across processes

    row = dict(pdb=pdb, n=int(t["n"]), fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
               prod_cloud=float(rec["rmsd_avg"]),
               pool_mean=float(rr_pool.mean()), pool_best=float(rr_pool.min()),
               set_mean={}, set_best={}, cells={})

    # ---- phase 1: ORACLE set-member statistics of the top-k (free)
    for k in K_GRID:
        idx = order[:k]
        row["set_mean"][str(k)] = float(rr_pool[idx].mean())
        row["set_best"][str(k)] = float(rr_pool[idx].min())

    # ---- phase 2: the actual operator -- coordinate average of a random m-subset of the top-k
    for k in K_GRID:
        top = order[:k]
        for m in M_GRID:
            if m > k:
                continue
            vals = []
            for s in range(N_SEED):
                S = top if m == k else top[rng.choice(k, m, replace=False)]
                Ps = P[np.ix_(S, S)]
                b = I.medoid(Ps)
                C = I.superpose_batch(W[S], W[S[b]]).mean(0)
                vals.append(I.ca_rmsd(C, nat))
            # deterministic deployable instantiation: every k/m-th member in SCORE RANK
            Sd = top[np.linspace(0, k - 1, m).round().astype(int)] if m < k else top
            Sd = np.unique(Sd)
            Pd = P[np.ix_(Sd, Sd)]
            Cd = I.superpose_batch(W[Sd], W[Sd[I.medoid(Pd)]]).mean(0)
            row["cells"]["k%d_m%d" % (k, m)] = dict(
                mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                stride=float(I.ca_rmsd(Cd, nat)), n_seed=len(vals))
    return row


def run_cloud():
    done = set()
    if os.path.exists(ROWS):
        with open(ROWS) as fh:
            for ln in fh:
                if ln.strip():
                    done.add(json.loads(ln)["pdb"])
    ts = I.targets()
    for t in ts:
        if t["pdb"] in done:
            continue
        row = one_target(t["pdb"], t)
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        print("%s  n=%2d  prod %.3f  setmean k75 %.3f k500 %.3f  cell k75m75 %.3f k500m75 %.3f"
              % (t["pdb"], t["n"], row["prod_cloud"], row["set_mean"]["75"], row["set_mean"]["500"],
                 row["cells"]["k75_m75"]["mean"], row["cells"]["k500_m75"]["mean"]), flush=True)
    print("cloud phase complete:", len(ts))


def analyse():
    R = [json.loads(ln) for ln in open(ROWS) if ln.strip()]
    assert len(R) == 126, len(R)
    pdbs = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])
    isf = np.array([r["fail18"] for r in R])
    pm = np.array([r["pool_mean"] for r in R])
    bp = np.array([r["pool_best"] for r in R])
    alt_pm = np.zeros(126, bool); alt_pm[np.argsort(-pm)[:18]] = True      # filter-independent tail 1
    alt_bp = np.zeros(126, bool); alt_bp[np.argsort(-bp)[:18]] = True      # filter-independent tail 2
    strata = dict(all=np.ones(126, bool), FAIL18=isf, other108=~isf,
                  worst18_poolmean=alt_pm, worst18_bestpool=alt_bp)

    out = dict(prereg="s30/PREREG_S30_F2.md", n=len(R),
               basis="POINT CLOUD; built chain is phase 3", K_GRID=list(K_GRID), M_GRID=list(M_GRID))

    # ---------- phase 1: ORACLE set-member statistics vs filter width
    sm = {k: np.array([r["set_mean"][str(k)] for r in R]) for k in K_GRID}
    sb = {k: np.array([r["set_best"][str(k)] for r in R]) for k in K_GRID}
    out["phase1_set_mean_ORACLE"] = {nm: {str(k): float(sm[k][s].mean()) for k in K_GRID}
                                     for nm, s in strata.items()}
    out["phase1_set_best_ORACLE"] = {nm: {str(k): float(sb[k][s].mean()) for k in K_GRID}
                                     for nm, s in strata.items()}
    # F2a clauses, measured against k = 500 (the no-filter end)
    ben = {nm: {str(k): float((sm[500][s] - sm[k][s]).mean()) for k in K_GRID} for nm, s in strata.items()}
    harm = {nm: {str(k): float((sb[k][s] - sb[500][s]).mean()) for k in K_GRID} for nm, s in strata.items()}
    out["phase1_setmean_benefit_vs_k500_ORACLE"] = ben
    out["phase1_setbest_harm_vs_k500_ORACLE"] = harm
    f2a = {}
    for k in K_GRID:
        if k <= 75:
            continue
        keep = ben["other108"][str(k)] / ben["other108"]["75"] if ben["other108"]["75"] else float("nan")
        shed = 1.0 - (harm["FAIL18"][str(k)] / harm["FAIL18"]["75"]) if harm["FAIL18"]["75"] else float("nan")
        f2a[str(k)] = dict(benefit_retained_108=float(keep), harm_shed_FAIL18=float(shed),
                           fires=bool(keep >= 0.80 and shed >= 0.50))
    out["F2a"] = dict(clauses=f2a, fires=any(v["fires"] for v in f2a.values()),
                      first_k_firing=next((k for k in sorted(f2a, key=int) if f2a[k]["fires"]), None))

    # ---------- phase 2: the real operator
    cell = {}
    for k in K_GRID:
        for m in M_GRID:
            if m <= k:
                cell[(k, m)] = np.array([r["cells"]["k%d_m%d" % (k, m)]["mean"] for r in R])
    out["phase2_cloud"] = {nm: {"k%d_m%d" % (k, m): float(v[s].mean()) for (k, m), v in cell.items()}
                           for nm, s in strata.items()}
    base = cell[(75, 75)]
    out["phase2_vs_production_m75"] = {}
    for k in K_GRID:
        c = ST.compare(cell[(k, 75)], base, folds=folds, names=pdbs,
                       label="k%d_m75_vs_k75m75" % k, seed_parts=("s30F2",))
        out["phase2_vs_production_m75"]["k%d" % k] = {
            q: c[q] for q in ("effect", "effect_over_mde", "mde", "ci95_fold", "folds_same_sign",
                              "n_folds", "n_better", "n_worse", "verdict", "mean_a", "mean_b")}
        for nm, s in strata.items():
            if nm == "all":
                continue
            cs = ST.compare(cell[(k, 75)][s], base[s], folds=folds[s],
                            names=[p for p, q in zip(pdbs, s) if q],
                            label="k%d_m75_%s" % (k, nm), seed_parts=("s30F2",))
            out["phase2_vs_production_m75"]["k%d" % k][nm] = {
                q: cs[q] for q in ("effect", "effect_over_mde", "ci95_fold", "folds_same_sign",
                                   "n_folds", "n_better", "n_worse", "verdict")}

    # ---------- F2c: is the k curve flat at fixed m = 75?
    mx = max(abs(out["phase2_vs_production_m75"]["k%d" % k]["effect_over_mde"]) for k in K_GRID)
    out["F2c"] = dict(max_abs_xMDE_over_k=float(mx), flat=bool(mx < 0.7),
                      note="flat => filter width behaves like averaging width; no-gate branch closes")

    # ---------- the averaging-width axis, for contrast with S29's flat m sweep
    out["averaging_width_at_fixed_k500"] = {
        "m%d" % m: float(cell[(500, m)].mean()) for m in M_GRID if m <= 500}
    out["averaging_width_at_fixed_k75"] = {
        "m%d" % m: float(cell[(75, m)].mean()) for m in M_GRID if m <= 75}

    out["provenance"] = ST.provenance(__file__)
    ST.save_atomic(OUT, out)

    print("\n=== PHASE 1  ORACLE set MEAN of the top-k (point cloud) ===")
    print("%-18s %s" % ("stratum", "".join("%8d" % k for k in K_GRID)))
    for nm in ("all", "other108", "FAIL18", "worst18_poolmean", "worst18_bestpool"):
        print("%-18s %s" % (nm, "".join("%8.3f" % out["phase1_set_mean_ORACLE"][nm][str(k)] for k in K_GRID)))
    print("\n=== PHASE 1  ORACLE set BEST of the top-k (point cloud) ===")
    for nm in ("all", "other108", "FAIL18", "worst18_poolmean", "worst18_bestpool"):
        print("%-18s %s" % (nm, "".join("%8.3f" % out["phase1_set_best_ORACLE"][nm][str(k)] for k in K_GRID)))
    print("\n=== F2a  benefit retained on the 108 / harm shed on FAIL18 (ORACLE) ===")
    for k in sorted(f2a, key=int):
        print("  k=%-4s benefit retained %6.1f%%   harm shed %6.1f%%   fires %s"
              % (k, 100 * f2a[k]["benefit_retained_108"], 100 * f2a[k]["harm_shed_FAIL18"], f2a[k]["fires"]))
    print("  F2a fires:", out["F2a"]["fires"], " first k:", out["F2a"]["first_k_firing"])
    print("\n=== PHASE 2  the real operator, POINT CLOUD, m = 75 fixed, k swept ===")
    print("%-6s %8s %9s %7s %-20s %6s %9s" % ("k", "mean", "vs prod", "xMDE", "foldCI", "folds", "W/L"))
    for k in K_GRID:
        c = out["phase2_vs_production_m75"]["k%d" % k]
        print("%-6d %8.4f %+9.4f %+7.2f [%+.4f,%+.4f] %4d/%d %5dW/%dL  %s"
              % (k, c["mean_a"], c["effect"], c["effect_over_mde"], c["ci95_fold"][0], c["ci95_fold"][1],
                 c["folds_same_sign"], c["n_folds"], c["n_better"], c["n_worse"], c["verdict"]))
    print("\n  per stratum (effect vs production, point cloud):")
    print("  %-6s %10s %10s %10s %10s" % ("k", "FAIL18", "other108", "w18_pmean", "w18_bpool"))
    for k in K_GRID:
        c = out["phase2_vs_production_m75"]["k%d" % k]
        print("  %-6d %+10.4f %+10.4f %+10.4f %+10.4f" % (
            k, c["FAIL18"]["effect"], c["other108"]["effect"],
            c["worst18_poolmean"]["effect"], c["worst18_bestpool"]["effect"]))
    print("\n  F2c max |xMDE| over k =", round(mx, 3), " flat:", out["F2c"]["flat"])
    print("\n  averaging width at k=500:", out["averaging_width_at_fixed_k500"])
    print("  averaging width at k=75 :", out["averaging_width_at_fixed_k75"])
    print("\nwrote", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["cloud", "analyse"])
    a = ap.parse_args()
    if a.phase == "cloud":
        run_cloud()
    else:
        analyse()
