#!/usr/bin/env python
"""s32/s32_Q2_shrinkctrl.py -- LANE Q: the control that decides whether hull projection is
anything other than shrinkage.

Q1-T2 consequence (2) says the convex readout beats emitting a structure estimate directly
only when the estimate's error exceeds the hull radius.  Before that can mean anything it needs
the control matched to the OPERATOR'S OWN SPACE (contract rule 6): projecting onto a bounded
convex set is a shrinkage, and this project has already priced shrinkage (S31 section 20.1's
shrink curve).  So each projection arm is compared against a shrinkage toward the pool mean
matched to **the projection's own displacement**, not to some other arm's.

Arms, per target, per noise level, per draw.  ALL **ORACLE / NOT DEPLOYABLE** (the estimate is
built from the native).  Basis: **CA point cloud**, a diagnostic.

    DIRECT    that = t + e                       emit that
    PROJ      x = P_conv{W}(that)                the convex readout
    SPAN      x = P_aff-span{W}(that)            projection onto the SPAN only (unbounded)
    SHRINK    x = Xbar + c (that - Xbar),  c chosen so ||x - Xbar|| == ||PROJ - Xbar||

If SHRINK matches PROJ, the hull contributes nothing beyond shrinkage and the observation is
not a mechanism.

    python s32/s32_Q2_shrinkctrl.py [--limit N] [--draws 8]
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

import core                                                            # noqa: E402
from core.pipeline import consensus_medoid                             # noqa: E402
from s8 import consensus2 as cc                                        # noqa: E402
from s32.s32_Q1_sufficiency import simplex_qp                          # noqa: E402

CACHE = os.path.join(ROOT, "s31", "results", "s31_C_cache")
OUT = os.path.join(HERE, "results")
NQ = 128
EPS = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--draws", type=int, default=8)
    args = ap.parse_args()

    geo = core.backend("geometry")
    names = sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz"))
    if args.limit:
        names = names[:args.limit]

    rows, t0 = [], time.time()
    for i, pdb in enumerate(names):
        rng = np.random.default_rng(abs(hash(("s32Q2", pdb))) % (2 ** 32))
        z = np.load(os.path.join(CACHE, "%s.npz" % pdb))
        order = np.asarray(z["order"], int)
        W = np.asarray(z["W"], float)
        nat = np.asarray(z["nat"], float)                              # ORACLE
        n, fold = int(z["n"]), int(z["fold"])
        Wo = W[order[:NQ]]
        Pt = geo.pairwise_ca_rmsd(Wo)
        ref = Wo[consensus_medoid(Pt)]
        U = cc.superpose_batch(Wo, ref).reshape(NQ, -1)
        t = cc.superpose_batch(nat[None], ref)[0].reshape(-1)          # ORACLE
        d = U.shape[1]
        scale = float(np.sqrt(n))
        G = U @ U.T
        Ub = U.mean(0)
        D = U - Ub
        _u, sv, vt = np.linalg.svd(D, full_matrices=False)
        rk = int((sv > sv[0] * 1e-10).sum())
        V = vt[:rk]

        per = []
        for eps in EPS:
            acc = {k: [] for k in ("direct", "proj", "span", "shrink", "cshr")}
            for _ in range(args.draws):
                e = rng.normal(size=d)
                e = eps * scale * e / np.linalg.norm(e)
                th = t + e
                w, _S, _k = simplex_qp(G, U @ th)
                xp = U.T @ w
                xs = Ub + V.T @ (V @ (th - Ub))                        # span only
                rp = np.linalg.norm(xp - Ub)
                rt = np.linalg.norm(th - Ub)
                c = rp / max(rt, 1e-30)                                # NORM-MATCHED to PROJ
                xk = Ub + c * (th - Ub)
                acc["direct"].append(float(geo.ca_rmsd(th.reshape(n, 3), nat)))
                acc["proj"].append(float(geo.ca_rmsd(xp.reshape(n, 3), nat)))
                acc["span"].append(float(geo.ca_rmsd(xs.reshape(n, 3), nat)))
                acc["shrink"].append(float(geo.ca_rmsd(xk.reshape(n, 3), nat)))
                acc["cshr"].append(float(c))
            per.append({"eps": float(eps),
                        **{k: float(np.mean(v)) for k, v in acc.items()},
                        **{k + "_sd": float(np.std(v, ddof=1)) for k, v in acc.items()}})
        rows.append({"pdb": pdb, "fold": fold, "n": n, "rank_aff": rk, "per": per})
        if i % 20 == 0:
            print("  %3d/%d %s  %.1fs" % (i + 1, len(names), pdb, time.time() - t0), flush=True)

    os.makedirs(OUT, exist_ok=True)
    tmp = os.path.join(OUT, "s32_Q2_shrinkctrl_rows.jsonl.%d.tmp" % os.getpid())
    with open(tmp, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    os.replace(tmp, os.path.join(OUT, "s32_Q2_shrinkctrl_rows.jsonl"))

    sys.path.insert(0, os.path.join(ROOT, "s24"))
    import stats_lib as S
    nm = [r["pdb"] for r in rows]
    fo = [r["fold"] for r in rows]

    table, cmps = [], []
    for k, eps in enumerate(EPS):
        g = lambda key: np.array([r["per"][k][key] for r in rows])     # noqa: E731
        table.append({"eps": eps, "direct": float(g("direct").mean()),
                      "proj": float(g("proj").mean()), "span": float(g("span").mean()),
                      "shrink": float(g("shrink").mean()),
                      "proj_se": float(g("proj").std(ddof=1) / np.sqrt(len(rows))),
                      "shrink_c_mean": float(g("cshr").mean()),
                      "draw_sd_mean": float(np.mean([r["per"][k]["proj_sd"] for r in rows]))})
        cmps.append(S.compare(g("proj"), g("direct"), folds=fo, names=nm,
                              label="PROJ-DIRECT_eps%.1f" % eps))
        cmps.append(S.compare(g("proj"), g("shrink"), folds=fo, names=nm,
                              label="PROJ-SHRINK_eps%.1f" % eps))

    summary = {"prereg_commit": "a8f9d6a7", "n": len(rows), "draws": args.draws,
               "label": "ORACLE / NOT DEPLOYABLE", "basis": "CA point cloud (diagnostic)",
               "table": table,
               "compare": {c["label"]: {kk: c[kk] for kk in
                                        ("effect", "se", "mde", "effect_over_mde",
                                         "n_better", "n_worse", "median_effect", "verdict")
                                        if kk in c} for c in cmps},
               "elapsed_s": time.time() - t0}
    tmp = os.path.join(OUT, "s32_Q2_shrinkctrl.json.%d.tmp" % os.getpid())
    with open(tmp, "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    os.replace(tmp, os.path.join(OUT, "s32_Q2_shrinkctrl.json"))

    print("%-6s %8s %8s %8s %8s %7s" % ("eps", "DIRECT", "PROJ", "SPAN", "SHRINK", "c"))
    for r in table:
        print("%-6.1f %8.4f %8.4f %8.4f %8.4f %7.3f"
              % (r["eps"], r["direct"], r["proj"], r["span"], r["shrink"], r["shrink_c_mean"]))
    for c in cmps:
        print("%-22s eff %+.4f SE %.4f MDE %.4f %.2fx %dW/%dL"
              % (c["label"], c["effect"], c["se"], c["mde"], c["effect_over_mde"],
                 c["n_better"], c["n_worse"]))


if __name__ == "__main__":
    main()
