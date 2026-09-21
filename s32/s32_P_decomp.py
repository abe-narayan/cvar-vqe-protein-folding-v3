#!/usr/bin/env python
"""s32/s32_P_decomp.py -- S32 LANE P, H-P1: the six-way pool decomposition, ONE STAGE AT A TIME.

PREREG s32/PREREG_S32_P.md @ 33dfe0d3.  Charter section 10 asks for POOL QUALITY vs POOL DIVERSITY
vs POOL COMMON-MODE BIAS vs WITHIN-POOL RANKING vs READOUT vs RECONSTRUCTION, **measured
independently**.  A cumulative ladder cannot do that: each rung is charged with the residual of
every earlier one, so the last rung always looks largest.

THE CELL NOBODY HAS MEASURED, and the reason this file exists: what a PERFECT RANKER buys
**through the operator production actually deploys**.  The S29 ladder's `best single member` rung
(2.3055 chain) changes the ranker AND the readout at once -- it reports the best member, which
production never emits.  Here the ranker is idealised and the readout is left exactly as shipped:

    ORACLE top-m by `rr`  ->  medoid-superposed uniform average  ->  (cloud)

ORACLE / NOT DEPLOYABLE for every idealised arm.  Basis: CA POINT CLOUD.  The chain price for
each arm is quoted from the S29 ladder where that arm exists there and is NOT invented where it
does not (contract rule 16: the projection price is a property of the object).

    python s32/s32_P_decomp.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                                        # noqa: E402
from s24 import stats_lib as ST                                       # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s32_P_decomp.json")
ROWS = os.path.join(RESULTS, "s32_P_decomp_rows.jsonl")
K, M, SEED = 500, 75, 32_0_8831
MS = [1, 3, 5, 10, 25, 50, 75]
NCOMP = [0]


def one(t):
    pdb = t["pdb"]
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    pool = np.asarray(u["order"], int)[:K]
    W = np.asarray(u["W"], float)[pool]
    rr = np.asarray(u["rr"], float)[pool]
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(float)
    sorder = np.argsort(DIS, kind="stable")
    oorder = np.argsort(rr, kind="stable")                     # ORACLE ranking of the same pool

    def emit(S):
        Ws = W[S]
        if len(S) == 1:
            return dict(cloud=float(I.ca_rmsd(Ws[0], nat)), setmean=float(rr[S].mean()),
                        setbest=float(rr[S].min()), wBw=0.0, f_common=1.0)
        P = I.pairwise_rmsd(Ws)
        b = I.medoid(P)
        Sup = I.superpose_batch(Ws, Ws[b])
        C = Sup.mean(0)
        X = Sup.reshape(len(S), -1)
        tf = I.superpose_batch(nat[None], Ws[b])[0].reshape(-1)
        e = X - tf[None, :]
        d = X - X.mean(0)[None, :]
        ebar2 = float(((X.mean(0) - tf) ** 2).sum())
        return dict(cloud=float(I.ca_rmsd(C, nat)), setmean=float(rr[S].mean()),
                    setbest=float(rr[S].min()), wBw=float(2.0 * (d ** 2).sum(1).mean()),
                    f_common=float(ebar2 / max((e ** 2).sum(1).mean(), 1e-12)))

    row = dict(pdb=pdb, n=int(t["n"]), fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
               pool_best=float(rr.min()), pool_mean=float(rr.mean()))
    for m in MS:
        for tag, o in (("SCORE", sorder), ("ORACLE", oorder)):
            for k, v in emit(o[:m]).items():
                row["%s_%s_m%d" % (k, tag, m)] = v
    # the shipped 75-set's own diversity secondaries, for the POOL DIVERSITY cell
    Ws = W[sorder[:M]]
    row["n_distinct_PROD"] = int(len({w.tobytes() for w in Ws}))
    P = I.pairwise_rmsd(Ws)
    row["pair_rmsd_mean_PROD"] = float(P[np.triu_indices(M, 1)].mean())
    row["pair_rmsd_sd_PROD"] = float(P[np.triu_indices(M, 1)].std())
    return row


def main():
    done = {}
    if os.path.exists(ROWS):
        for ln in open(ROWS):
            if ln.strip():
                r = json.loads(ln); done[r["pdb"]] = r
    ts = I.targets()
    R = []
    with open(ROWS, "a") as fh:
        for t in ts:
            if t["pdb"] in done:
                R.append(done[t["pdb"]]); continue
            r = one(t); R.append(r); fh.write(json.dumps(r) + "\n"); fh.flush()
    R.sort(key=lambda r: r["pdb"])
    assert len(R) == 126 and len({r["pdb"] for r in R}) == 126, (len(R), len({r["pdb"] for r in R}))
    names = [r["pdb"] for r in R]; folds = np.array([r["fold"] for r in R])
    g = lambda k: np.array([r[k] for r in R], float)                       # noqa: E731

    def cmp2(a, b, lab):
        NCOMP[0] += 1
        o = ST.compare(a, b, folds=folds, names=names, label=lab, seed_parts=("s32Pdec", str(SEED)))
        return {k: o[k] for k in ("label", "mean_a", "mean_b", "effect", "median_effect", "se",
                                  "mde", "effect_over_mde", "ci95_fold", "folds_same_sign",
                                  "n_better", "n_worse", "verdict") if k in o}

    prod = g("cloud_SCORE_m75")
    curves = {}
    for tag in ("SCORE", "ORACLE"):
        curves[tag] = {m: dict(cloud=float(g("cloud_%s_m%d" % (tag, m)).mean()),
                               setmean=float(g("setmean_%s_m%d" % (tag, m)).mean()),
                               setbest=float(g("setbest_%s_m%d" % (tag, m)).mean()),
                               wBw=float(g("wBw_%s_m%d" % (tag, m)).mean()),
                               f_common=float(g("f_common_%s_m%d" % (tag, m)).mean()))
                       for m in MS}

    # ---------------------------------------------------------------- the one-at-a-time cells
    cells = {}
    cells["WITHIN_POOL_RANKING"] = dict(
        what="ORACLE top-75 by rr, THEN the unchanged medoid-superposed uniform average",
        value=float(g("cloud_ORACLE_m75").mean()),
        vs_production=cmp2(g("cloud_ORACLE_m75"), prod,
                           "ORACLE top-75 through the DEPLOYED operator - production (CLOUD)"),
        note="this is the cell the S29 ladder does not have: it idealises the RANKER and leaves "
             "the READOUT shipped.  The ladder's `best single member` changes both at once.")
    best_m = min(MS, key=lambda m: curves["ORACLE"][m]["cloud"])
    cells["WITHIN_POOL_RANKING_best_m_ORACLE"] = dict(
        what="ORACLE top-m through the deployed average, m chosen ORACLE-globally",
        m=best_m, value=curves["ORACLE"][best_m]["cloud"],
        vs_production=cmp2(g("cloud_ORACLE_m%d" % best_m), prod,
                           "ORACLE top-%d through the DEPLOYED operator - production (CLOUD)"
                           % best_m))
    cells["POOL_DIVERSITY"] = dict(
        n_distinct_of_75=float(g("n_distinct_PROD").mean()),
        mean_pairwise_rmsd=float(g("pair_rmsd_mean_PROD").mean()),
        sd_pairwise_rmsd=float(g("pair_rmsd_sd_PROD").mean()),
        wBw=float(g("wBw_SCORE_m75").mean()))
    cells["POOL_COMMON_MODE"] = dict(
        f_common_on_the_shipped_75=float(g("f_common_SCORE_m75").mean()),
        f_common_on_the_ORACLE_75=float(g("f_common_ORACLE_m75").mean()),
        note="f is descriptive only and is NOT a screen (S24 correction); |ebar|^2 = n*RMSD^2 "
             "exactly, so the numerator is not independent evidence")
    cells["POOL_QUALITY"] = dict(
        pool_best_K500=float(g("pool_best").mean()), pool_mean_K500=float(g("pool_mean").mean()),
        note="the universe-wide best (1.3134) and the matched-random controls live in "
             "s32_P_ladder.json; generation headroom on the BEST axis is 1.7108 -> 1.3134")

    out = dict(prereg="s32/PREREG_S32_P.md @ 33dfe0d3", hypothesis="H-P1",
               basis="CA POINT CLOUD; ORACLE / NOT DEPLOYABLE for every ORACLE_* arm",
               n=126, ms=MS, production_cloud=float(prod.mean()),
               curve_SCORE_prefix=curves["SCORE"], curve_ORACLE_prefix=curves["ORACLE"],
               cells=cells,
               operator_law=None, multiplicity_emitted=int(NCOMP[0]))

    # the operator's own law, fitted across BOTH prefix families and all m (S18 says a
    # score-based gate breaks it, so it is fitted here rather than quoted)
    A, y = [], []
    for tag in ("SCORE", "ORACLE"):
        for m in MS:
            if m == 1:
                continue
            A.append(np.column_stack([np.ones(126), g("setmean_%s_m%d" % (tag, m)),
                                      g("setbest_%s_m%d" % (tag, m))]))
            y.append(g("cloud_%s_m%d" % (tag, m)))
    A = np.vstack(A); y = np.concatenate(y)
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    out["operator_law"] = dict(
        basis="CLOUD, both prefix families, m in 3..75",
        intercept=float(beta[0]), coef_set_mean=float(beta[1]), coef_set_best=float(beta[2]),
        ratio=float(beta[1] / beta[2]) if beta[2] else None,
        R2=float(1 - ((A @ beta - y) ** 2).sum() / ((y - y.mean()) ** 2).sum()),
        registered_prediction="coef_set_mean >= 5x coef_set_best")
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
