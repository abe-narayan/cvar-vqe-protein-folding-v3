#!/usr/bin/env python
"""s32/s32_P_spread.py -- S32 LANE P: the untested middle between "no quality" and "no spread".

Coordinator's arm.  The readout objective is `<w,a> - 1/2 w'Bw`: it rewards LOW quality cost and
HIGH spread.  Two extremes are measured -- quality-blind dispersion maximisation (+0.1436 A, 1.22x,
WORSE, S31) and production (a score-ordered prefix with no spread term at all).  The middle:

    use the score only as a FLOOR -- drop the worst fraction, the one thing it demonstrably does
    well -- then choose the 75 to MAXIMISE w'Bw among the survivors, uniform weights.

ONE DEGENERACY, reported because it removes a control the request asked for: "floor-only with a
score-ordered prefix among survivors" is IDENTICAL TO PRODUCTION for every floor that leaves at
least 75 survivors, because the survivors ARE the best-scoring prefix.  The non-degenerate floor
control is a RANDOM 75 from the same survivor set, which is run here with its own draw
distribution (rule 10).

SCREEN FIRST.  All arms are screened on the CA POINT CLOUD, which is free; only an arm that clears
1 MDE on the cloud is worth 17 projections per target.  Stated in advance so a cloud number is
never quoted as an endpoint.

    python s32/s32_P_spread.py screen
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

from s12 import instrument as I                                        # noqa: E402
from s24 import stats_lib as ST                                       # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s32_P_spread_rows.jsonl")
OUT = os.path.join(RESULTS, "s32_P_spread.json")
K, M, SEED, NDRAW = 500, 75, 32_0_4477, 5
FLOORS = [0.25, 0.50, 0.60, 0.75]          # fixed a priori; NOT swept finely (rule 11)
NCOMP = [0]


def greedy_maxdisp(B, m):
    """Uniform-weight maximiser of w'Bw over subsets of size m: the standard 2-approx greedy
    (seed with the most distant pair, then add the point maximising the sum of B to the set)."""
    nb = len(B)
    i, j = np.unravel_index(int(np.argmax(B)), B.shape)
    S = [int(i), int(j)]
    tot = B[:, S].sum(1).copy()
    inS = np.zeros(nb, bool); inS[S] = True
    while len(S) < m:
        tot_m = np.where(inS, -np.inf, tot)
        k = int(np.argmax(tot_m))
        S.append(k); inS[k] = True
        tot = tot + B[:, k]
    return np.array(S, int)


def one(t):
    pdb = t["pdb"]; n = int(t["n"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    pool = np.asarray(u["order"], int)[:K]
    W = np.asarray(u["W"], float)[pool]
    rr = np.asarray(u["rr"], float)[pool]
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(float)
    order = np.argsort(DIS, kind="stable")

    def emit(S):
        Ws = W[S]
        P = I.pairwise_rmsd(Ws)
        b = I.medoid(P)
        Sup = I.superpose_batch(Ws, Ws[b])
        X = Sup.reshape(len(S), -1)
        d = X - X.mean(0)[None, :]
        return dict(cloud=float(I.ca_rmsd(Sup.mean(0), nat)),
                    setmean=float(rr[S].mean()), setbest=float(rr[S].min()),
                    wBw=float(2.0 * (d ** 2).sum(1).mean()))

    row = dict(pdb=pdb, n=n, fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18))
    for k, v in emit(order[:M]).items():
        row[k + "_PROD"] = v
    for f in FLOORS:
        keep = order[:max(M, int(round(K * (1.0 - f))))]
        Wk = W[keep]
        Pk = I.pairwise_rmsd(Wk)
        bk = I.medoid(Pk)
        Xk = I.superpose_batch(Wk, Wk[bk]).reshape(len(keep), -1)
        sq = (Xk ** 2).sum(1)
        B = np.maximum(sq[:, None] + sq[None, :] - 2.0 * (Xk @ Xk.T), 0.0)
        np.fill_diagonal(B, 0.0)
        S = keep[greedy_maxdisp(B, M)]
        tag = "F%02d" % int(round(100 * f))
        row["nkeep_" + tag] = int(len(keep))
        for k, v in emit(S).items():
            row["%s_SPREAD_%s" % (k, tag)] = v
        for dd in range(NDRAW):
            g = np.random.default_rng(SEED + 31 * dd + (hash(pdb + tag) & 0xFFFF))
            Sr = keep[g.choice(len(keep), M, replace=False)]
            for k, v in emit(Sr).items():
                row["%s_RANDFLOOR_%s_%d" % (k, tag, dd)] = v
    return row


def screen():
    done = set()
    if os.path.exists(ROWS):
        for ln in open(ROWS):
            if ln.strip():
                done.add(json.loads(ln)["pdb"])
    ts = I.targets(); t0 = time.time()
    for k, t in enumerate(ts):
        if t["pdb"] in done:
            continue
        r = one(t)
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        if (k + 1) % 10 == 0:
            print("[%3d/126 %5.0fs] %s cloud PROD %.3f SPREAD_F50 %.3f RANDFLOOR_F50 %.3f"
                  % (k + 1, time.time() - t0, r["pdb"], r["cloud_PROD"],
                     r["cloud_SPREAD_F50"],
                     np.mean([r["cloud_RANDFLOOR_F50_%d" % d] for d in range(NDRAW)])), flush=True)
    analyse()


def analyse():
    R = [json.loads(ln) for ln in open(ROWS) if ln.strip()]
    R.sort(key=lambda r: r["pdb"])
    if len(R) != 126:
        raise RuntimeError("have %d rows, not 126" % len(R))
    names = [r["pdb"] for r in R]; folds = np.array([r["fold"] for r in R])
    g = lambda k: np.array([r[k] for r in R], float)                       # noqa: E731

    def cmp2(a, b, lab):
        NCOMP[0] += 1
        o = ST.compare(a, b, folds=folds, names=names, label=lab, seed_parts=("s32Psp", str(SEED)))
        return {k: o[k] for k in ("label", "mean_a", "mean_b", "effect", "median_effect", "se",
                                  "mde", "effect_over_mde", "ci95_fold", "folds_same_sign",
                                  "n_better", "n_worse", "verdict") if k in o}

    prod = g("cloud_PROD")
    out = dict(prereg="s32/PREREG_S32_P.md @ 33dfe0d3 (coordinator arm, EXPLORATORY)",
               basis="CA POINT CLOUD -- a SCREEN, not an endpoint.  Nothing here is a chain result.",
               n=126, floors=FLOORS, ndraw=NDRAW, seed=SEED,
               degeneracy="floor + score-ordered prefix == PRODUCTION exactly for every floor "
                          "leaving >= 75 survivors; that control is not run because it is the "
                          "incumbent under another name",
               PROD=dict(cloud=float(prod.mean()), setmean=float(g("setmean_PROD").mean()),
                         setbest=float(g("setbest_PROD").mean()), wBw=float(g("wBw_PROD").mean())),
               arms={})
    for f in FLOORS:
        tag = "F%02d" % int(round(100 * f))
        sp = g("cloud_SPREAD_" + tag)
        rf = np.stack([g("cloud_RANDFLOOR_%s_%d" % (tag, d)) for d in range(NDRAW)], 1)
        out["arms"][tag] = dict(
            floor=f, n_survivors=float(g("nkeep_" + tag).mean()),
            SPREAD=dict(cloud=float(sp.mean()), setmean=float(g("setmean_SPREAD_" + tag).mean()),
                        setbest=float(g("setbest_SPREAD_" + tag).mean()),
                        wBw=float(g("wBw_SPREAD_" + tag).mean()),
                        vs_PROD=cmp2(sp, prod, "SPREAD %s - PROD (CLOUD SCREEN)" % tag)),
            RANDFLOOR=dict(cloud_draw_mean=float(rf.mean(0).mean()),
                           cloud_draw_sd=float(rf.mean(0).std(ddof=1)),
                           setmean=float(np.mean([g("setmean_RANDFLOOR_%s_%d" % (tag, d)).mean()
                                                  for d in range(NDRAW)])),
                           wBw=float(np.mean([g("wBw_RANDFLOOR_%s_%d" % (tag, d)).mean()
                                              for d in range(NDRAW)])),
                           per_draw_vs_PROD=[cmp2(rf[:, d], prod,
                                                  "RANDFLOOR %s draw %d - PROD (CLOUD SCREEN)"
                                                  % (tag, d)) for d in range(NDRAW)]),
            SPREAD_vs_RANDFLOOR=cmp2(sp, rf.mean(1),
                                     "SPREAD %s - RANDFLOOR mean-of-draws (CLOUD SCREEN) -- "
                                     "THE CONTROL THAT DECIDES WHETHER SPREAD DOES ANY WORK" % tag))
    #: leave-fold-out floor choice, on TRAINING-FOLD natives only
    sp = {f: g("cloud_SPREAD_F%02d" % int(round(100 * f))) for f in FLOORS}
    lfo = np.empty(126)
    for fo in sorted(set(folds.tolist())):
        tr = folds != fo
        best = min(FLOORS, key=lambda f: sp[f][tr].mean())
        lfo[folds == fo] = sp[best][folds == fo]
        out.setdefault("LFO_floor_by_fold", {})[int(fo)] = best
    out["SPREAD_LFO_vs_PROD"] = cmp2(lfo, prod, "SPREAD, floor chosen leave-fold-out - PROD (CLOUD SCREEN)")
    out["GATE"] = ("PROJECT" if out["SPREAD_LFO_vs_PROD"]["effect_over_mde"] <= -1.0 else
                   "DO NOT PROJECT -- no arm clears 1 MDE on the free screen")
    out["multiplicity_emitted"] = int(NCOMP[0])
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["screen", "analyse"])
    a = ap.parse_args()
    (screen if a.phase == "screen" else analyse)()
