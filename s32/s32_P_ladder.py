#!/usr/bin/env python
"""s32/s32_P_ladder.py -- S32 LANE P, H-P2a: the filter ladder with best-of-K priced at EVERY rung.

PREREG: s32/PREREG_S32_P.md @ 33dfe0d3.  Contract rule 9 (a per-target minimum over K variants is
mostly an order statistic) and rule 10 (a random control needs its own draw distribution).

ORACLE / NOT DEPLOYABLE everywhere: every `best member` reads `rr`, the native CA-RMSD label.

THE RUNG NOBODY HAS PRICED.  The S29 ladder starts at the K=500 pool.  That pool is itself a
FILTER: `core/pipeline.py:698-700` keeps the 500 highest BLOSUM62 sums out of a window universe of
7,016 to 39,410 members.  The 500 -> 128 cut is a DIFFERENT filter -- the distogram Bayes-risk
score (`core/pipeline.py:755-760`) -- and the 128 window exists only when `quantum=True`, which
production is not.  So this file measures four rungs and names the filter at each.

Each rung reports, at MATCHED K, both quantities the pipeline can actually spend:
    BEST  member of the retained set   (needs a perfect ranker; the S29 ladder's axis)
    MEAN  over the retained set        (what the uniform-average terminal operator consumes)

    python s32/s32_P_ladder.py
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
OUT = os.path.join(RESULTS, "s32_P_ladder.json")
SEED = 32_0_5501
NDRAW = 200                     # rule 10: a distribution, never a draw
NCOMP = [0]


def main():
    ts = I.targets()
    folds = np.array([t["fold"] for t in ts]); names = [t["pdb"] for t in ts]
    rung = {}
    per = {k: [] for k in
           ("UNI_best", "UNI_mean", "nw",
            "B500_best", "B500_mean", "R500u_best", "R500u_mean", "R500u_best_sd", "R500u_mean_sd",
            "S128_best", "S128_mean", "R128_best", "R128_mean", "R128_best_sd", "R128_mean_sd",
            "S75_best", "S75_mean", "R75of128_best", "R75of128_mean",
            "R75of500_best", "R75of500_mean", "R75of500_best_sd", "R75of500_mean_sd",
            "B75_best", "B75_mean",
            "rank_of_pool_best_in_score", "pool_best_survives_128", "pool_best_survives_75")}

    for t in ts:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        rr = np.asarray(u["rr"], float)
        order = np.asarray(u["order"], int)              # BLOSUM ordering of the whole universe
        nw = len(rr)
        pool = order[:500]
        rp = rr[pool]
        DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(float)
        so = np.argsort(DIS, kind="stable")

        rng = np.random.default_rng(SEED + abs(hash(pdb)) % 100003)

        def draws(src, k, n=NDRAW):
            idx = np.array([rng.choice(len(src), k, replace=False) for _ in range(n)])
            v = src[idx]
            return v.min(1), v.mean(1)

        b, m = draws(rr, 500)
        per["R500u_best"].append(b.mean()); per["R500u_best_sd"].append(b.std(ddof=1))
        per["R500u_mean"].append(m.mean()); per["R500u_mean_sd"].append(m.std(ddof=1))
        b, m = draws(rp, 128)
        per["R128_best"].append(b.mean()); per["R128_best_sd"].append(b.std(ddof=1))
        per["R128_mean"].append(m.mean()); per["R128_mean_sd"].append(m.std(ddof=1))
        b, m = draws(rp[so[:128]], 75)
        per["R75of128_best"].append(b.mean()); per["R75of128_mean"].append(m.mean())
        b, m = draws(rp, 75)
        per["R75of500_best"].append(b.mean()); per["R75of500_best_sd"].append(b.std(ddof=1))
        per["R75of500_mean"].append(m.mean()); per["R75of500_mean_sd"].append(m.std(ddof=1))

        per["nw"].append(nw)
        per["UNI_best"].append(rr.min()); per["UNI_mean"].append(rr.mean())
        per["B500_best"].append(rp.min()); per["B500_mean"].append(rp.mean())
        per["S128_best"].append(rp[so[:128]].min()); per["S128_mean"].append(rp[so[:128]].mean())
        per["S75_best"].append(rp[so[:75]].min()); per["S75_mean"].append(rp[so[:75]].mean())
        #: the 75 highest BLOSUM sums -- the retrieval key used as the gate as well
        per["B75_best"].append(rp[:75].min()); per["B75_mean"].append(rp[:75].mean())
        kbest = int(np.argmin(rp))
        per["rank_of_pool_best_in_score"].append(int(np.nonzero(so == kbest)[0][0]) + 1)
        per["pool_best_survives_128"].append(float(per["rank_of_pool_best_in_score"][-1] <= 128))
        per["pool_best_survives_75"].append(float(per["rank_of_pool_best_in_score"][-1] <= 75))

    A = {k: np.asarray(v, float) for k, v in per.items()}

    def cmp2(a, b, lab):
        NCOMP[0] += 1
        o = ST.compare(a, b, folds=folds, names=names, label=lab, seed_parts=("s32Plad", str(SEED)))
        return {k: o[k] for k in ("label", "mean_a", "mean_b", "effect", "median_effect", "se",
                                  "mde", "effect_over_mde", "ci95_fold", "folds_same_sign",
                                  "n_better", "n_worse", "verdict") if k in o}

    def rungrec(name, filt, parent, child, rnd_b, rnd_m, rnd_b_sd=None, rnd_m_sd=None):
        pb, cb = A[parent + "_best"], A[child + "_best"]
        pm, cm = A[parent + "_mean"], A[child + "_mean"]
        r = dict(
            filter=filt,
            BEST=dict(parent=float(pb.mean()), child=float(cb.mean()),
                      rung_loss=float((cb - pb).mean()),
                      matched_random_at_child_K=float(rnd_b.mean()),
                      order_statistic_part=float((rnd_b - pb).mean()),
                      attributable_to_the_filter=float((cb - rnd_b).mean()),
                      filter_vs_random=cmp2(cb, rnd_b, "%s BEST: filter - matched random" % name),
                      draw_sd_within_target=float(rnd_b_sd.mean()) if rnd_b_sd is not None else None),
            MEAN=dict(parent=float(pm.mean()), child=float(cm.mean()),
                      rung_gain=float((cm - pm).mean()),
                      matched_random_at_child_K=float(rnd_m.mean()),
                      order_statistic_part=float((rnd_m - pm).mean()),
                      attributable_to_the_filter=float((cm - rnd_m).mean()),
                      filter_vs_random=cmp2(cm, rnd_m, "%s MEAN: filter - matched random" % name),
                      draw_sd_within_target=float(rnd_m_sd.mean()) if rnd_m_sd is not None else None))
        rung[name] = r
        return r

    rungrec("R1_universe_to_500", "BLOSUM62 sum (core/pipeline.py:698-700)",
            "UNI", "B500", A["R500u_best"], A["R500u_mean"], A["R500u_best_sd"], A["R500u_mean_sd"])
    rungrec("R2_500_to_128", "distogram Bayes risk (core/pipeline.py:755-760) -- QUANTUM ONLY",
            "B500", "S128", A["R128_best"], A["R128_mean"], A["R128_best_sd"], A["R128_mean_sd"])
    rungrec("R3_128_to_75", "distogram Bayes risk prefix", "S128", "S75",
            A["R75of128_best"], A["R75of128_mean"])
    rungrec("R2b_500_to_75_DEPLOYABLE", "distogram Bayes risk -- THE FILTER PRODUCTION ACTUALLY RUNS",
            "B500", "S75", A["R75of500_best"], A["R75of500_mean"],
            A["R75of500_best_sd"], A["R75of500_mean_sd"])

    out = dict(
        prereg="s32/PREREG_S32_P.md @ 33dfe0d3", hypothesis="H-P2a", n=126, ndraw=NDRAW, seed=SEED,
        basis="CA POINT CLOUD quality labels rr; ORACLE / NOT DEPLOYABLE (every `best` reads rr)",
        universe=dict(n_windows_min=int(A["nw"].min()), n_windows_median=float(np.median(A["nw"])),
                      n_windows_max=int(A["nw"].max()), best=float(A["UNI_best"].mean()),
                      mean=float(A["UNI_mean"].mean())),
        ladder_best=dict(universe=float(A["UNI_best"].mean()), blosum500=float(A["B500_best"].mean()),
                         score128=float(A["S128_best"].mean()), score75=float(A["S75_best"].mean()),
                         blosum75=float(A["B75_best"].mean())),
        ladder_mean=dict(universe=float(A["UNI_mean"].mean()), blosum500=float(A["B500_mean"].mean()),
                         score128=float(A["S128_mean"].mean()), score75=float(A["S75_mean"].mean()),
                         blosum75=float(A["B75_mean"].mean())),
        rungs=rung,
        survival=dict(rank_of_pool_best_in_score_mean=float(A["rank_of_pool_best_in_score"].mean()),
                      rank_median=float(np.median(A["rank_of_pool_best_in_score"])),
                      survives_128=float(A["pool_best_survives_128"].mean()),
                      survives_75=float(A["pool_best_survives_75"].mean()),
                      uninformative_null_survives_128=128.0 / 500.0,
                      uninformative_null_survives_75=75.0 / 500.0),
        blosum75_vs_score75=dict(
            BEST=cmp2(A["B75_best"], A["S75_best"], "BLOSUM top-75 - score top-75, BEST"),
            MEAN=cmp2(A["B75_mean"], A["S75_mean"], "BLOSUM top-75 - score top-75, MEAN")),
        multiplicity_emitted=int(NCOMP[0]))

    # ------------------------------------------------------------------ CONCENTRATION, and the
    # stratum that makes the BEST-axis headline unquotable as a broad effect.  The coordinator
    # retracted "the score is WORSE than random" on exactly this ground: `n_better` counts
    # targets where the SCORE wins, so 72W/54L means the score wins on 72 of 126 while losing
    # more on the ones it loses.  `median-vs-mean-is-the-free-warning` says print both, plus a
    # uniform-effect null, as ONE verdict.
    f18 = np.array([p in I.FAIL18 for p in names], bool)
    conc = {}
    for nm, child, rnd in (("R2_500_to_128_BEST", A["S128_best"], A["R128_best"]),
                           ("R2b_500_to_75_BEST", A["S75_best"], A["R75of500_best"]),
                           ("R1_universe_to_500_BEST", A["B500_best"], A["R500u_best"])):
        d_ = child - rnd
        c = ST.concentration(d_, seed_parts=("s32Plad", nm))
        conc[nm] = dict(mean=float(d_.mean()), median=float(np.median(d_)),
                        W_score_better=int((d_ < 0).sum()), L_score_worse=int((d_ > 0).sum()),
                        concentration=c,
                        FAIL18_CIRCULAR_mean=float(d_[f18].mean()),
                        other108_mean=float(d_[~f18].mean()),
                        other108_over_mde=float(d_[~f18].mean() /
                                                (2.8016 * d_[~f18].std(ddof=1) / np.sqrt((~f18).sum()))),
                        share_of_total_from_FAIL18=float(d_[f18].sum() / d_.sum()) if d_.sum() else None,
                        note="FAIL18 is defined by no pool member within 1.5 A of the pool best "
                             "surviving the top-75 filter, so a filter-vs-random contrast on it "
                             "is NEAR-CIRCULAR; the other-108 row is the quotable one")
    out["concentration_and_strata"] = conc
    out["per_target"] = {k: [float(v) for v in A[k]] for k in
                         ("UNI_best", "B500_best", "S128_best", "S75_best", "R500u_best",
                          "R128_best", "R75of500_best", "UNI_mean", "B500_mean", "S128_mean",
                          "S75_mean", "R500u_mean", "R128_mean", "R75of500_mean")}
    out["names"] = names
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
