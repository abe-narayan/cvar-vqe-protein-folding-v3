#!/usr/bin/env python
"""s31/s31_V_orderstat.py -- LANE V (adversary), V1: the MATCHED-K control for the
set-matched readout ladder (STATE NOTE 2(b) / S31-L6 section 1).

FALSIFIER, REGISTERED IN s31/AUDIT_V.md BEFORE THIS RAN
-------------------------------------------------------
STATE NOTE 2(b) reads the fixed-top-128 ladder as

    argmin over 128            7.0 bits             2.1458
    2-of-128, UNIFORM weights  12.99 support bits   2.0700   "-0.076 vs argmin, ZERO weight bits"

and the coordinator's handed-back prediction is that the 0.076 A is ERROR CANCELLATION.

But `2-of-128` is a per-target minimum over C(128,2) = 8128 ORACLE supports, while `argmin over
128` is a per-target minimum over 128.  That is a 64x larger oracle search, and project memory
`grid-oracles-are-order-statistics` plus contract rule 11 say a per-target minimum over K variants
is mostly best-of-K.  Lane F has just measured (S31-L11) that on this very object the per-target
minimum is MONOTONE AND UNSATURATED in K.

V1 REGISTERED BAR, written before the number exists:
  "the 0.076 A is an order statistic, not error cancellation" FIRES if the per-target minimum over
  **128 RANDOM pairs** (matched K, matched operator, same top-128, uniform weights) fails to beat
  the per-target minimum over the 128 singletons by at least 0.7 x MDE -- i.e. if you need the
  full 8128-pair search to get the 0.076.

Everything here is CA POINT CLOUD and ORACLE / NOT DEPLOYABLE.  No chain claim is made.
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

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s31_V_orderstat_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_V_orderstat.json")
SEED = 31099
K128 = 128
NDRAW = 8            # contract rule 10: a draw control needs its own distribution; MEAN, never max
KGRID = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8128]
NCOMP = [0]


def top128(pdb):
    """The deployed top-128 of the 500-member pool under the shipped score, plus the native."""
    u = I.load_univ(pdb)
    dg = I.distogram(pdb)
    n = int(u["n"])
    pool = np.asarray(u["order"], int)[:500]
    W = np.asarray(u["W"], float)[pool]
    nat = np.asarray(u["nat_ca"], float)
    ii, jj = I.pair_index(n, 2)
    D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    order = np.argsort(I.shipped_score(dg, D), kind="stable")
    return W[order[:K128]], nat, n


def run():
    done = set()
    if os.path.exists(ROWS):
        for ln in open(ROWS):
            if ln.strip():
                done.add(json.loads(ln)["pdb"])
    iu, ju = np.triu_indices(K128, 1)               # 8128 pairs, fixed order
    for t in I.targets():
        pdb = t["pdb"]
        if pdb in done:
            continue
        Wt, nat, n = top128(pdb)
        singles = I.kabsch_rmsd_batch(Wt, nat)      # 128 -- min is best1_top128 (CLOUD)
        # all 8128 UNIFORM 2-member averages.  ca_rmsd is Kabsch, so the anchor choice inside the
        # pair is irrelevant: two anchorings differ by a rigid transform.
        pairs = np.empty(len(iu))
        for a in range(K128 - 1):
            sel = ju[iu == a]
            if not len(sel):
                continue
            moved = I.superpose_batch(Wt[sel], Wt[a])
            avg = 0.5 * (moved + Wt[a][None])
            pairs[iu == a] = I.kabsch_rmsd_batch(avg, nat)
        row = dict(pdb=pdb, n=n, fold=int(t["fold"]),
                   singles_min=float(singles.min()), singles_mean=float(singles.mean()),
                   pairs_min=float(pairs.min()), pairs_mean=float(pairs.mean()),
                   singles=singles.tolist(), pairs_min_by_anchor=None)
        # order-statistic growth: min over K RANDOM pairs, and min over K RANDOM singles
        rng = np.random.default_rng(SEED + (abs(int.from_bytes(pdb.encode(), "big")) % 100003))
        gp, gs = {}, {}
        for K in KGRID:
            if K <= len(pairs):
                gp[str(K)] = float(np.mean([pairs[rng.choice(len(pairs), K, replace=False)].min()
                                            for _ in range(NDRAW)]))
            if K <= K128:
                gs[str(K)] = float(np.mean([singles[rng.choice(K128, K, replace=False)].min()
                                            for _ in range(NDRAW)]))
        row["growth_pairs"] = gp
        row["growth_singles"] = gs
        # the matched-K arm, kept per draw so the draw distribution is reportable (rule 10)
        row["pairs_min_K128_per_draw"] = [
            float(pairs[rng.choice(len(pairs), K128, replace=False)].min()) for _ in range(NDRAW)]
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        print("%s n=%2d  best1_128 %.3f  best-pair(8128) %.3f  best-pair(K=128 rand) %.3f"
              % (pdb, n, row["singles_min"], row["pairs_min"],
                 float(np.mean(row["pairs_min_K128_per_draw"]))), flush=True)
    print("V1 rows complete")


def cmp2(a, b, folds, names, label):
    NCOMP[0] += 1
    return ST.compare(np.asarray(a, float), np.asarray(b, float), folds=folds, names=names,
                      label=label, seed_parts=("s31V", str(SEED)))


def brief(o):
    return dict(label=o["label"], mean_a=o["mean_a"], mean_b=o["mean_b"], effect=o["effect"],
                median_effect=o["median_effect"], se=o["se"], mde=o["mde"],
                effect_over_mde=o["effect_over_mde"], ci95_fold=o["ci95_fold"],
                folds_same_sign=o["folds_same_sign"], W=o["n_better"], L=o["n_worse"], n=o["n"])


def gate(o):
    r = abs(o["effect_over_mde"]); ci = o["ci95_fold"]
    ex = bool(ci is not None and (ci[0] > 0 or ci[1] < 0))
    if r >= 1.0 and ex:
        return "MEASURED"
    if r >= 1.0:
        return "NOT MEASURED (>=1.0x MDE, fold CI spans zero)"
    if r >= 0.7:
        return "NOT MEASURED (0.7-1.0x MDE)"
    return "NOT A RESULT (<0.7x MDE)"


def analyse():
    R = [json.loads(ln) for ln in open(ROWS) if ln.strip()]
    R.sort(key=lambda r: r["pdb"])
    assert len(R) == 126, "INCOMPLETE: %d/126 -- contract rule 15, check the job" % len(R)
    names = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])
    s_min = np.array([r["singles_min"] for r in R])           # = best1_top128, CLOUD
    p_min = np.array([r["pairs_min"] for r in R])             # exhaustive 2-of-128 uniform
    p_k128 = np.array([np.mean(r["pairs_min_K128_per_draw"]) for r in R])
    per_draw = np.array([r["pairs_min_K128_per_draw"] for r in R])   # 126 x NDRAW

    out = dict(seed=SEED, n=len(R), basis="CA POINT CLOUD",
               ORACLE="every arm here is ORACLE / NOT DEPLOYABLE -- each is a per-target minimum "
                      "taken against the native",
               prereg="s31/AUDIT_V.md V1, registered before this script was run",
               reproduction=dict(
                   best1_top128_recomputed=float(s_min.mean()),
                   best1_top128_s31C_ladder=2.1457974841561676,
                   note="lane C's ORACLE_argmin128; agreement gates everything below"))
    out["reproduction"]["max_abs_dev_vs_C"] = abs(float(s_min.mean()) - 2.1457974841561676)

    o_ex = cmp2(p_min, s_min, folds, names,
                "V1.exhaustive 2-of-128 uniform (K=8128) - argmin over 128 (K=128), CLOUD")
    o_mk = cmp2(p_k128, s_min, folds, names,
                "V1.MATCHED-K 2-of-128 uniform (K=128 random pairs) - argmin over 128, CLOUD")
    out["V1a_exhaustive_pairs_vs_argmin"] = dict(cmp=brief(o_ex), gate=gate(o_ex))
    out["V1b_MATCHED_K_pairs_vs_argmin"] = dict(
        cmp=brief(o_mk), gate=gate(o_mk), n_draws=NDRAW,
        per_draw_mean=[float(per_draw[:, d].mean()) for d in range(per_draw.shape[1])],
        per_draw_sd=float(np.std([per_draw[:, d].mean() for d in range(per_draw.shape[1])], ddof=1)),
        note="MEAN over draws, never the maximum (contract rule 10)")

    g_ex = float(p_min.mean() - s_min.mean())
    g_mk = float(p_k128.mean() - s_min.mean())
    out["V1_REGISTERED_BAR"] = dict(
        bar="'the 0.076 A is an order statistic, not error cancellation' FIRES if the matched-K "
            "pair family fails to beat the 128 singletons by >= 0.7 x MDE",
        exhaustive_gain_vs_argmin=g_ex, matched_K_gain_vs_argmin=g_mk,
        share_of_exhaustive_gain_surviving_at_matched_K=(g_mk / g_ex) if g_ex else float("nan"),
        BAR_FIRES=bool(g_mk > -0.7 * o_mk["mde"]))

    gp = {K: float(np.mean([r["growth_pairs"][str(K)] for r in R]))
          for K in KGRID if str(KGRID[0]) in R[0]["growth_pairs"] and str(K) in R[0]["growth_pairs"]}
    gs = {K: float(np.mean([r["growth_singles"][str(K)] for r in R]))
          for K in KGRID if str(K) in R[0]["growth_singles"]}
    out["V1c_order_statistic_growth"] = dict(
        pairs_min_over_K=gp, singles_min_over_K=gs,
        note="mean over 126 targets of the per-target minimum over K randomly chosen variants; "
             "NDRAW draws per cell, mean over draws")

    # --- the bit-budget axis, against S30's OWN reference curve (s30_Q_sparse.json :: argmin_ref)
    qs = os.path.join(ROOT, "s30", "results", "s30_Q_sparse.json")
    if os.path.exists(qs):
        Q = json.load(open(qs))
        ref = {int(k): v["mean"] for k, v in Q["argmin_ref"].items()}
        arms = {k: v for k, v in Q["arms"].items() if k.startswith("T128_")}
        rows = []
        for k, v in sorted(arms.items()):
            wb = v.get("w_bits")
            wb = float(wb) if wb not in (None, "Infinity") else float("inf")
            sup = float(v.get("sup_bits", float("nan")))
            tot = sup + wb
            rows.append(dict(arm=k, mean=v["mean"], sup_bits=sup, w_bits=wb, total_bits=tot))
        out["V1d_bit_budget_vs_S30_own_reference"] = dict(
            argmin_ref=ref, T128_arms=rows,
            reading="S30-L11's claim is 'a plain argmin dominates at EVERY BIT BUDGET'. On S30's "
                    "own argmin_ref curve the plain argmin reaches 1.9383 at 8 bits and 1.7108 at "
                    "9 bits, against T128_s2_unif's 2.0700 at 12.99 support bits. The claim is "
                    "NOT inverted by the set-matched ladder; what lane C showed is that S30 "
                    "SUPPORTED a correct claim with a set-mismatched pair.")
    out["multiplicity"] = dict(family="s31_V_orderstat", comparisons_emitted=int(NCOMP[0]))
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "V1d_bit_budget_vs_S30_own_reference"},
                     indent=1, default=float)[:6000])
    print("\nwrote", OUT)


if __name__ == "__main__":
    ph = sys.argv[1] if len(sys.argv) > 1 else "run"
    if ph == "run":
        run()
    else:
        analyse()
