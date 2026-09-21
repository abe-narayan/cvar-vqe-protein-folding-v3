#!/usr/bin/env python
"""s32/s32_P_rand.py -- S32 LANE P, H-P2b: the DEPLOYABLE random-gate arms on the BUILT CHAIN.

Coordinator's request (S32-L?, lane V's order-statistic control): lane V showed on the CA cloud
that the score's top-128 holds a WORSE best member than a random 128 of the same K=500 pool
(+0.1872, 1.06x MDE).  That is a statement about the best MEMBER.  Production does not take the
best member -- it averages 75.  This job asks the only question that decides anything:

    does replacing the score's gate with a RANDOM gate of the same size lower the BUILT CHAIN?

A random gate needs no native information, so a win here is immediately DEPLOYABLE and needs no
in-band skill, which is the wall every other route has hit.

ARMS, per target, ALL PROJECTED IN THE SAME PROCESS AS PROD (contract rule 3):

    PROD     the score's top-75 of the K=500 pool               (the deployable incumbent)
    R128_d   a random 128 of the 500, then the SCORE's top-75 within it   d = 0..NDRAW-1
    R75_d    a random 75 of the 500, the score never consulted           d = 0..NDRAW-1

Every arm also reports its SET MEAN and SET BEST (ORACLE, diagnostic) and the uniform readout's
dispersion w'Bw = 2*mean_k|d_k|^2, so a change can be attributed to quality or to diversity, plus
the virtual-bond mean/sd the contract requires of any arm from its first row (rule 15).

    python s32/s32_P_rand.py run [--shard i --of N]
    python s32/s32_P_rand.py analyse
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from zlib import crc32

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
ROWPAT = os.path.join(RESULTS, "s32_P_rand_rows.s%dof%d.jsonl")
OUT = os.path.join(RESULTS, "s32_P_rand.json")
SEED = 32_0_2119
NDRAW = 8
K, K128, M = 500, 128, 75
NCOMP = [0]


def stats_of(Wsel, rr_sel, nat):
    """The uniform-average readout on a selected set: cloud, set quality, dispersion."""
    P = I.pairwise_rmsd(Wsel)
    b = I.medoid(P)
    Sup = I.superpose_batch(Wsel, Wsel[b])
    C = Sup.mean(0)
    X = Sup.reshape(len(Wsel), -1)
    d = X - X.mean(0)[None, :]
    return C, dict(cloud=float(I.ca_rmsd(C, nat)),
                   set_mean=float(np.mean(rr_sel)), set_best=float(np.min(rr_sel)),
                   set_median=float(np.median(rr_sel)),
                   mean_d2=float((d ** 2).sum(1).mean()),
                   wBw=float(2.0 * (d ** 2).sum(1).mean()),
                   n_distinct=int(len({w.tobytes() for w in Wsel})))


def one(t):
    pdb = t["pdb"]
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    pool = np.asarray(u["order"], int)[:K]
    W = np.asarray(u["W"], float)[pool]
    rr = np.asarray(u["rr"], float)[pool]
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(float)
    order = np.argsort(DIS, kind="stable")

    sets = {"PROD": order[:M]}
    for dd in range(NDRAW):
        g128 = np.random.default_rng(SEED + 101 * dd + crc32((pdb + "|128").encode()))
        sel = g128.choice(K, K128, replace=False)
        # the score's top-75 WITHIN the random 128 -- the pipeline's own filter, new gate
        sets["R128_%d" % dd] = sel[np.argsort(DIS[sel], kind="stable")[:M]]
        g75 = np.random.default_rng(SEED + 907 * dd + crc32((pdb + "|75").encode()))
        sets["R75_%d" % dd] = g75.choice(K, M, replace=False)

    row = dict(pdb=pdb, n=int(t["n"]), fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
               pool_best=float(rr.min()), pool_mean=float(rr.mean()))
    for nm, S in sets.items():
        S = np.asarray(S, int)
        C, st = stats_of(W[S], rr[S], nat)
        pr = I.project(C, t["seq"], t["fold"])
        ca = np.asarray(pr["ca"], float)
        bl = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        row["cloud_" + nm] = st["cloud"]
        row["chain_" + nm] = float(I.ca_rmsd(ca, nat))
        row["setmean_" + nm] = st["set_mean"]
        row["setbest_" + nm] = st["set_best"]
        row["setmed_" + nm] = st["set_median"]
        row["wBw_" + nm] = st["wBw"]
        row["ndist_" + nm] = st["n_distinct"]
        row["bond_" + nm] = [float(bl.mean()), float(bl.std())]
    return row


def run(shard, nshard):
    path = ROWPAT % (shard, nshard)
    done = set()
    if os.path.exists(path):
        for ln in open(path):
            if ln.strip():
                done.add(json.loads(ln)["pdb"])
    ts = [t for k, t in enumerate(I.targets()) if k % nshard == shard]
    t0 = time.time()
    for k, t in enumerate(ts):
        if t["pdb"] in done:
            continue
        r = one(t)
        with open(path, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        rc = np.mean([r["chain_R128_%d" % d] for d in range(NDRAW)])
        r7 = np.mean([r["chain_R75_%d" % d] for d in range(NDRAW)])
        print("[%3d/%3d %6.0fs] %s  chain PROD %6.3f  R128 %6.3f  R75 %6.3f | setmean PROD %6.3f "
              "R75 %6.3f" % (k + 1, len(ts), time.time() - t0, r["pdb"], r["chain_PROD"], rc, r7,
                             r["setmean_PROD"], np.mean([r["setmean_R75_%d" % d] for d in range(NDRAW)])),
              flush=True)
    print("shard %d/%d complete" % (shard, nshard))


def _rows():
    R = []
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_P_rand_rows.s*of*.jsonl"))):
        for ln in open(f):
            if ln.strip():
                R.append(json.loads(ln))
    seen, out = set(), []
    for r in sorted(R, key=lambda x: x["pdb"]):
        if r["pdb"] not in seen:
            seen.add(r["pdb"]); out.append(r)
    return out


def analyse():
    R = _rows()
    if len(R) != 126:
        raise RuntimeError("have %d distinct rows, not 126 -- the run is not finished" % len(R))
    names = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])
    g = lambda k: np.array([r[k] for r in R], float)                       # noqa: E731

    def cmp2(a, b, lab):
        NCOMP[0] += 1
        o = ST.compare(a, b, folds=folds, names=names, label=lab, seed_parts=("s32P", str(SEED)))
        return {k: o[k] for k in ("label", "n", "mean_a", "mean_b", "median_a", "median_b",
                                  "effect", "median_effect", "se", "mde", "effect_over_mde",
                                  "ci95_fold", "folds_same_sign", "n_better", "n_worse",
                                  "verdict") if k in o}

    # ------------------------------------------------------------ rule 3: print the check
    ref = {}
    for f in glob.glob(os.path.join(ROOT, "s29", "results", "s29_O_chain_rows*.jsonl")):
        for ln in open(f):
            if ln.strip():
                d = json.loads(ln)
                if d.get("item") == "prod":
                    ref.setdefault(d["pdb"], d["rmsd_chain"])
    prod = g("chain_PROD")
    have = np.array([p in ref for p in names])
    dref = np.array([prod[i] - ref[names[i]] for i in range(len(names)) if have[i]])
    repro = dict(n_matched=int(have.sum()), mean_abs=float(np.abs(dref).mean()),
                 p90_abs=float(np.percentile(np.abs(dref), 90)), max_abs=float(np.abs(dref).max()),
                 mean_signed=float(dref.mean()), my_PROD_mean=float(prod.mean()),
                 s29_PROD_mean=float(np.mean([ref[n] for n in names if n in ref])),
                 A2_floor="0.0134 mean / 0.0329 p90 / 0.2285 max (contract rule 3)",
                 within_A2_floor=bool(np.abs(dref).mean() <= 0.0134))

    out = dict(prereg="s32/PREREG_S32_P.md @ 33dfe0d3", hypothesis="H-P2b (coordinator arm)",
               basis="BUILT CHAIN, every arm projected in the same process as PROD for that target",
               n=126, ndraw=NDRAW, seed=SEED, PROD_reproduction=repro)

    for fam, lab in (("R128", "random 128 of 500, then the SCORE's top-75 within"),
                     ("R75", "random 75 of 500, the score never consulted")):
        ch = np.stack([g("chain_%s_%d" % (fam, d)) for d in range(NDRAW)], 1)
        cl = np.stack([g("cloud_%s_%d" % (fam, d)) for d in range(NDRAW)], 1)
        sm = np.stack([g("setmean_%s_%d" % (fam, d)) for d in range(NDRAW)], 1)
        sb = np.stack([g("setbest_%s_%d" % (fam, d)) for d in range(NDRAW)], 1)
        wb = np.stack([g("wBw_%s_%d" % (fam, d)) for d in range(NDRAW)], 1)
        per = [cmp2(ch[:, d], prod, "%s draw %d - PROD (BUILT CHAIN)" % (fam, d))
               for d in range(NDRAW)]
        out[fam] = dict(
            what=lab,
            chain_draw_means=[float(v) for v in ch.mean(0)],
            chain_draw_mean=float(ch.mean(0).mean()), chain_draw_sd=float(ch.mean(0).std(ddof=1)),
            cloud_draw_mean=float(cl.mean(0).mean()), cloud_draw_sd=float(cl.mean(0).std(ddof=1)),
            setmean_draw_mean=float(sm.mean(0).mean()), setbest_draw_mean=float(sb.mean(0).mean()),
            wBw_draw_mean=float(wb.mean(0).mean()),
            per_draw_vs_PROD=per,
            effect_draw_mean=float(np.mean([o["effect"] for o in per])),
            effect_draw_sd=float(np.std([o["effect"] for o in per], ddof=1)),
            ratio_draw_mean=float(np.mean([o["effect_over_mde"] for o in per])),
            n_draws_better=int(sum(o["effect"] < 0 for o in per)),
            AVG_OF_DRAWS_vs_PROD=cmp2(ch.mean(1), prod,
                                      "%s per-target mean over %d draws - PROD (BUILT CHAIN)"
                                      % (fam, NDRAW)))

    out["PROD"] = dict(chain=float(prod.mean()), cloud=float(g("cloud_PROD").mean()),
                       setmean=float(g("setmean_PROD").mean()),
                       setbest=float(g("setbest_PROD").mean()),
                       wBw=float(g("wBw_PROD").mean()),
                       ndist=float(g("ndist_PROD").mean()),
                       pool_best=float(g("pool_best").mean()),
                       pool_mean=float(g("pool_mean").mean()))

    # ------------------------------------------------------------ the operator's own law
    # d_out ~ set_mean + set_best, re-measured across ALL arms (S18: the law BREAKS under a
    # score-based gate, so it is fitted here and not quoted).
    A, y, tag = [], [], []
    for nm in ["PROD"] + ["R128_%d" % d for d in range(NDRAW)] + ["R75_%d" % d for d in range(NDRAW)]:
        A.append(np.column_stack([np.ones(126), g("setmean_" + nm), g("setbest_" + nm)]))
        y.append(g("chain_" + nm)); tag += [nm.split("_")[0]] * 126
    A = np.vstack(A); y = np.concatenate(y); tag = np.array(tag)
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    r2 = float(1 - ((A @ beta - y) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    law = dict(all_arms=dict(intercept=float(beta[0]), coef_set_mean=float(beta[1]),
                             coef_set_best=float(beta[2]), R2=r2,
                             ratio_mean_over_best=float(beta[1] / beta[2]) if beta[2] else None))
    for fam in ("PROD", "R128", "R75"):
        m = tag == fam
        b2, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
        law[fam] = dict(intercept=float(b2[0]), coef_set_mean=float(b2[1]),
                        coef_set_best=float(b2[2]),
                        R2=float(1 - ((A[m] @ b2 - y[m]) ** 2).sum() / ((y[m] - y[m].mean()) ** 2).sum()))
    out["operator_law"] = law
    out["geometry"] = {nm: dict(bond_mean=float(np.mean([r["bond_" + nm][0] for r in R])),
                                bond_sd=float(np.mean([r["bond_" + nm][1] for r in R])))
                       for nm in ("PROD", "R128_0", "R75_0")}
    out["multiplicity_emitted"] = int(NCOMP[0])
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["run", "analyse"])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    a = ap.parse_args()
    if a.phase == "run":
        run(a.shard, a.of)
    else:
        analyse()
