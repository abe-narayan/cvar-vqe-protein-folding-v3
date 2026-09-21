#!/usr/bin/env python
"""s32/s32_P_rankchain.py -- S32 LANE P: H-P1's biggest cell, carried to the BUILT CHAIN.

`s32_P_decomp.py` measures what a PERFECT RANKER buys through the operator production actually
deploys, and finds -1.0853 A at m=75 and -1.5670 A at the ORACLE-global m=5 -- on the CA CLOUD.
Lane P's own rules say a cloud result does not count until it reaches the chain, and contract
rule 16 says the projection price is a property of the OBJECT (~0 for a valid chain, ~+0.16 for
a dense average), so a 75-member average and a 5-member average are NOT in the same regime and
the price cannot be assumed equal.  All three arms are therefore projected IN THE SAME PROCESS.

    PROD      the score's top-75              (the deployable incumbent)
    ORA75     the ORACLE top-75 by rr         ORACLE / NOT DEPLOYABLE
    ORA5      the ORACLE top-5 by rr          ORACLE / NOT DEPLOYABLE (m fixed at 5 from the
                                              cloud analysis, not re-chosen on the chain)

    python s32/s32_P_rankchain.py run [--shard i --of N]
    python s32/s32_P_rankchain.py analyse
"""
from __future__ import annotations

import argparse
import glob
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
ROWPAT = os.path.join(RESULTS, "s32_P_rankchain_rows.s%dof%d.jsonl")
OUT = os.path.join(RESULTS, "s32_P_rankchain.json")
K, M, MSTAR, SEED = 500, 75, 5, 32_0_7703
NCOMP = [0]


def one(t):
    pdb = t["pdb"]
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    pool = np.asarray(u["order"], int)[:K]
    W = np.asarray(u["W"], float)[pool]
    rr = np.asarray(u["rr"], float)[pool]
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(float)
    sets = {"PROD": np.argsort(DIS, kind="stable")[:M],
            "ORA75": np.argsort(rr, kind="stable")[:M],
            "ORA5": np.argsort(rr, kind="stable")[:MSTAR]}
    row = dict(pdb=pdb, n=int(t["n"]), fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18))
    for nm, S in sets.items():
        Ws = W[S]
        P = I.pairwise_rmsd(Ws)
        b = I.medoid(P)
        C = I.superpose_batch(Ws, Ws[b]).mean(0)
        pr = I.project(C, t["seq"], t["fold"])
        ca = np.asarray(pr["ca"], float)
        bl = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        row["cloud_" + nm] = float(I.ca_rmsd(C, nat))
        row["chain_" + nm] = float(I.ca_rmsd(ca, nat))
        row["setmean_" + nm] = float(rr[S].mean())
        row["setbest_" + nm] = float(rr[S].min())
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
        print("[%3d/%3d %6.0fs] %s chain PROD %6.3f ORA75 %6.3f ORA5 %6.3f | cloud %6.3f %6.3f %6.3f"
              % (k + 1, len(ts), time.time() - t0, r["pdb"], r["chain_PROD"], r["chain_ORA75"],
                 r["chain_ORA5"], r["cloud_PROD"], r["cloud_ORA75"], r["cloud_ORA5"]), flush=True)
    print("shard %d/%d complete" % (shard, nshard))


def analyse():
    raw = []
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_P_rankchain_rows.s*of*.jsonl"))):
        for ln in open(f):
            if ln.strip():
                raw.append(json.loads(ln))
    by, conflict = {}, []
    for r in raw:
        if r["pdb"] in by and json.dumps(by[r["pdb"]], sort_keys=True) != json.dumps(r, sort_keys=True):
            conflict.append(r["pdb"])
        by[r["pdb"]] = r
    print("rows: %d raw, %d distinct, %d conflicting" % (len(raw), len(by), len(conflict)))
    if conflict:
        raise RuntimeError("duplicate pdbs with differing content: %s" % conflict[:5])
    R = [by[p] for p in sorted(by)]
    if len(R) != 126:
        raise RuntimeError("have %d rows, not 126" % len(R))
    names = [r["pdb"] for r in R]; folds = np.array([r["fold"] for r in R])
    f18 = np.array([r["fail18"] for r in R], bool)
    g = lambda k: np.array([r[k] for r in R], float)                       # noqa: E731

    def cmp2(a, b, lab, mask=None):
        NCOMP[0] += 1
        m = np.ones(len(folds), bool) if mask is None else np.asarray(mask, bool)
        o = ST.compare(np.asarray(a, float), np.asarray(b, float), folds=folds[m],
                       names=[n for n, k in zip(names, m) if k], label=lab,
                       seed_parts=("s32Prc", str(SEED)))
        return {k: o[k] for k in ("label", "mean_a", "mean_b", "median_a", "median_b", "effect",
                                  "median_effect", "se", "mde", "effect_over_mde", "ci95_fold",
                                  "folds_same_sign", "n_better", "n_worse", "verdict") if k in o}

    # rule 3: cross-job chain comparison is legal only if the shared PROD rows are checked
    ref = {}
    for f in glob.glob(os.path.join(ROOT, "s29", "results", "s29_O_chain_rows*.jsonl")):
        for ln in open(f):
            if ln.strip():
                d = json.loads(ln)
                if d.get("item") == "prod":
                    ref.setdefault(d["pdb"], d["rmsd_chain"])
    dref = np.array([g("chain_PROD")[i] - ref[names[i]] for i in range(126) if names[i] in ref])

    out = dict(prereg="s32/PREREG_S32_P.md @ 33dfe0d3", hypothesis="H-P1 on the BUILT CHAIN",
               basis="BUILT CHAIN, all three arms projected in the same process per target",
               ORACLE="ORA75 and ORA5 are ORACLE / NOT DEPLOYABLE (they rank by rr)",
               n=126, m_star_from_cloud=MSTAR,
               PROD_reproduction=dict(n_matched=int(len(dref)), mean_abs=float(np.abs(dref).mean()),
                                      p90_abs=float(np.percentile(np.abs(dref), 90)),
                                      max_abs=float(np.abs(dref).max()),
                                      my_mean=float(g("chain_PROD").mean()),
                                      s29_mean=float(np.mean([ref[n] for n in names if n in ref])),
                                      A2_floor="0.0134 mean / 0.0329 p90 / 0.2285 max"),
               means={k: dict(cloud=float(g("cloud_" + k).mean()), chain=float(g("chain_" + k).mean()),
                              projection_price=float((g("chain_" + k) - g("cloud_" + k)).mean()),
                              setmean=float(g("setmean_" + k).mean()),
                              setbest=float(g("setbest_" + k).mean()),
                              bond_mean=float(np.mean([r["bond_" + k][0] for r in R])),
                              bond_sd=float(np.mean([r["bond_" + k][1] for r in R])))
                      for k in ("PROD", "ORA75", "ORA5")},
               ORA75_vs_PROD=cmp2(g("chain_ORA75"), g("chain_PROD"),
                                  "ORACLE top-75 through the DEPLOYED average - PROD (BUILT CHAIN)"),
               ORA5_vs_PROD=cmp2(g("chain_ORA5"), g("chain_PROD"),
                                 "ORACLE top-5 through the DEPLOYED average - PROD (BUILT CHAIN)"),
               ORA5_vs_ORA75=cmp2(g("chain_ORA5"), g("chain_ORA75"),
                                  "ORACLE top-5 - ORACLE top-75 (BUILT CHAIN)"),
               strata=dict(
                   ORA75_vs_PROD_other108=cmp2(g("chain_ORA75")[~f18], g("chain_PROD")[~f18],
                                               "ORA75 - PROD, other 108 (BUILT CHAIN)", mask=~f18),
                   ORA75_vs_PROD_FAIL18_CIRCULAR=cmp2(g("chain_ORA75")[f18], g("chain_PROD")[f18],
                                                      "ORA75 - PROD, FAIL18 (BUILT CHAIN) "
                                                      "[CIRCULAR STRATUM]", mask=f18)),
               multiplicity_emitted=int(NCOMP[0]))
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
