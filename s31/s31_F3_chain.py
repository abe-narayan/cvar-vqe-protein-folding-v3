#!/usr/bin/env python
"""s31/s31_F3_chain.py -- S31 LANE F, F3 on the BUILT CHAIN, all arms in ONE JOB.

The 146% matched-random-family result of S31-L11 is a CA-point-cloud comparison.  This carries it
to the endpoint.  Contract: lane D showed the lam=0 multi-start argmin is decided at a 1e-7 spread
between branches 1e-1 apart (amplification ~1e13), so cross-job chain comparisons are unsafe --
**every arm here is projected in the same process from the same stored clouds.**

    M75      prefix average of the top-75 of the top-128  (the production operator)
    PREFIX   prefix average at the CLOUD-argmin m         (= bestm128, as S29 built it)
    RANDOM   the random-subset average at the CLOUD-argmin variant, matched family, draw 0
    RANDOM1..3  the same for draws 1..3                    (contract rule 10: a draw distribution)

Every arm's variant is selected on the CLOUD, exactly as `s29_O_ladder.py:542` selects `m_best`,
so the two oracles are matched in what they are allowed to see.  ORACLE / NOT DEPLOYABLE, all.

    python s31/s31_F3_chain.py run
    python s31/s31_F3_chain.py analyse
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

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
ROWS = os.path.join(RESULTS, "s31_F3_chain_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_F3_chain.json")
SEED, K128, NDRAW = 31007, 128, 2
NCOMP = [0]
PREFJSON = os.path.join(RESULTS, "s31_F3_prefix.json")


def m_settings():
    """The ORACLE-global m and the leave-fold-out m, read from the cloud analysis that derived
    them (`s31_F3_prefix.json` -> F3c), NOT re-derived here.  Both are S29-L30's own arms; this
    job exists only to put them on the BUILT CHAIN, because S29-L30's transfer arms are point
    cloud and the cloud->chain price on this rung is +0.1422."""
    d = json.load(open(PREFJSON))["F3c_global_m"]
    return int(d["m_star"]), {int(k): int(v) for k, v in d["m_by_fold"].items()}


def one(t, m_star, m_fold):
    pdb = t["pdb"]
    u = I.load_univ(pdb); dg = I.distogram(pdb); n = int(u["n"])
    pool = np.asarray(u["order"], int)[:500]
    W = np.asarray(u["W"], float)[pool]
    nat = np.asarray(u["nat_ca"], float)
    ii, jj = I.pair_index(n, 2)
    D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    order = np.argsort(I.shipped_score(dg, D), kind="stable")
    Wt = W[order[:K128]]
    Pt = I.pairwise_rmsd(Wt)

    def avg(S):
        b = I.medoid(Pt[np.ix_(S, S)])
        return I.superpose_batch(Wt[S], Wt[S[b]]).mean(0)

    pref = np.empty(K128)
    Cp = [None] * K128
    for m in range(1, K128 + 1):
        Cp[m - 1] = avg(np.arange(m))
        pref[m - 1] = I.ca_rmsd(Cp[m - 1], nat)
    mb = int(pref.argmin())

    mlfo = int(m_fold[int(t["fold"])])
    arms = {"M75": Cp[74], "M_GLOBAL": Cp[m_star - 1], "M_LFO": Cp[mlfo - 1],
            "PREFIX": Cp[mb]}
    rand_cloud = {}
    for d in range(NDRAW):
        rng = np.random.default_rng(SEED + 1000 * d + crc32(pdb.encode()))
        v = np.empty(K128); Cs = [None] * K128
        for m in range(1, K128 + 1):
            S = rng.choice(K128, m, replace=False) if m < K128 else np.arange(K128)
            Cs[m - 1] = avg(S)
            v[m - 1] = I.ca_rmsd(Cs[m - 1], nat)
        k = int(v.argmin())
        arms["RANDOM%d" % d] = Cs[k]
        rand_cloud["RANDOM%d" % d] = dict(cloud=float(v[k]), variant=k + 1)

    row = dict(pdb=pdb, n=n, fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
               m_best=mb + 1, m_star=int(m_star), m_lfo=mlfo,
               cloud_M75=float(pref[74]), cloud_M_GLOBAL=float(pref[m_star - 1]),
               cloud_M_LFO=float(pref[mlfo - 1]), cloud_PREFIX=float(pref[mb]),
               **{("cloud_" + k): v["cloud"] for k, v in rand_cloud.items()},
               **{("variant_" + k): v["variant"] for k, v in rand_cloud.items()})
    for nm, C in arms.items():
        pr = I.project(np.asarray(C, float), t["seq"], t["fold"])
        ca = np.asarray(pr["ca"], float)
        row["chain_" + nm] = float(I.ca_rmsd(ca, nat))
        d_ = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        row["bond_" + nm] = [float(d_.mean()), float(d_.std())]
    return row


def run():
    done = set()
    if os.path.exists(ROWS):
        for ln in open(ROWS):
            if ln.strip():
                done.add(json.loads(ln)["pdb"])
    ts = I.targets(); t0 = time.time()
    m_star, m_fold = m_settings()
    print("ORACLE global m* = %d ; leave-fold-out m = %s" % (m_star, m_fold), flush=True)
    for k, t in enumerate(ts):
        if t["pdb"] in done:
            continue
        r = one(t, m_star, m_fold)
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        print("[%3d/%3d %5.0fs] %s m*=%3d | cloud M75 %6.3f PREFIX %6.3f RAND0 %6.3f | "
              "chain M75 %6.3f PREFIX %6.3f RAND0 %6.3f"
              % (k + 1, len(ts), time.time() - t0, r["pdb"], r["m_best"],
                 r["cloud_M75"], r["cloud_PREFIX"], r["cloud_RANDOM0"],
                 r["chain_M75"], r["chain_PREFIX"], r["chain_RANDOM0"]), flush=True)
    print("run complete")


def analyse():
    R = [json.loads(ln) for ln in open(ROWS) if ln.strip()]
    R.sort(key=lambda r: r["pdb"])
    if len(R) != 126:
        raise RuntimeError("have %d rows, not 126 -- the run is not finished" % len(R))
    names = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])

    def c(k):
        return np.array([r[k] for r in R], float)

    def cmp2(a, b, lab):
        NCOMP[0] += 1
        o = ST.compare(a, b, folds=folds, names=names, label=lab,
                       seed_parts=("s31F3chain", str(SEED)))
        return dict(label=o["label"], mean_a=o["mean_a"], mean_b=o["mean_b"],
                    effect=o["effect"], median_effect=o["median_effect"], se=o["se"],
                    mde=o["mde"], effect_over_mde=o["effect_over_mde"],
                    ci95_fold=o["ci95_fold"], folds_same_sign=o["folds_same_sign"],
                    W=o["n_better"], L=o["n_worse"], n=o["n"])

    m75c, prefc = c("chain_M75"), c("chain_PREFIX")
    glob_c, lfo_c = c("chain_M_GLOBAL"), c("chain_M_LFO")
    rnd = np.stack([c("chain_RANDOM%d" % d) for d in range(NDRAW)], 1)
    out = dict(seed=SEED, n=126, basis="BUILT CHAIN, all arms projected in ONE job from the same "
                                       "stored clouds; variant selected on the CLOUD for every arm",
               ORACLE="ORACLE / NOT DEPLOYABLE -- every arm is a per-target minimum over K = 128",
               means=dict(M75=float(m75c.mean()), PREFIX=float(prefc.mean()),
                          M_GLOBAL=float(glob_c.mean()), M_LFO=float(lfo_c.mean()),
                          RANDOM_per_draw=[float(v) for v in rnd.mean(0)],
                          RANDOM_mean=float(rnd.mean(0).mean()),
                          RANDOM_sd_over_draws=float(rnd.mean(0).std(ddof=1))),
               PREFIX_vs_M75=cmp2(prefc, m75c, "F3chain.PREFIX(bestm128) - M75 (BUILT CHAIN)"),
               M_GLOBAL_vs_M75=cmp2(glob_c, m75c, "F3chain.ORACLE global m - M75 (BUILT CHAIN)"),
               M_LFO_vs_M75=cmp2(lfo_c, m75c, "F3chain.leave-fold-out m - M75 (BUILT CHAIN)"),
               RANDOM_vs_M75=[cmp2(rnd[:, d], m75c,
                                   "F3chain.RANDOM draw %d - M75 (BUILT CHAIN)" % d)
                              for d in range(NDRAW)],
               RANDOM_vs_PREFIX=[cmp2(rnd[:, d], prefc,
                                      "F3chain.RANDOM draw %d - PREFIX (BUILT CHAIN)" % d)
                                 for d in range(NDRAW)])
    gp = out["PREFIX_vs_M75"]["effect"]
    gr = float(np.mean([o["effect"] for o in out["RANDOM_vs_M75"]]))
    out["share_of_prefix_gain_BUILT_CHAIN"] = float(gr / gp) if gp else float("nan")
    out["BAR_FIRES_on_chain"] = bool(gp and (gr / gp) >= 0.80)
    out["cloud_cross_check"] = dict(
        M75=float(c("cloud_M75").mean()), PREFIX=float(c("cloud_PREFIX").mean()),
        M_GLOBAL=float(c("cloud_M_GLOBAL").mean()), M_LFO=float(c("cloud_M_LFO").mean()),
        cloud_to_chain_price=dict(M75=float((c("chain_M75")-c("cloud_M75")).mean()),
                                  PREFIX=float((prefc-c("cloud_PREFIX")).mean()),
                                  M_GLOBAL=float((glob_c-c("cloud_M_GLOBAL")).mean()),
                                  M_LFO=float((lfo_c-c("cloud_M_LFO")).mean())),
        RANDOM_mean=float(np.mean([c("cloud_RANDOM%d" % d).mean() for d in range(NDRAW)])),
        note="must agree with s31_F3_prefix.json's cloud figures (2.7605 prefix) up to the "
             "reseeding of the random family")
    out["geometry"] = {k: dict(bond_mean=float(np.mean([r["bond_" + k][0] for r in R])),
                               bond_sd=float(np.mean([r["bond_" + k][1] for r in R])))
                       for k in ("M75", "PREFIX", "RANDOM0")}
    out["multiplicity_emitted"] = int(NCOMP[0])
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["run", "analyse"])
    a = ap.parse_args()
    (run if a.phase == "run" else analyse)()
