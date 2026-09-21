#!/usr/bin/env python
"""s31/s31_C_widen.py -- lane C, question C3 of `s31/PREREG_S31_C.md`.

Discharges the circularity gate on S30-L11's incidental finding that widening the quantum register
128 -> 512 is worth -1.9004 A on FAIL18.

EVERY NUMBER THIS FILE PRODUCES IS **ORACLE / NOT DEPLOYABLE**.  `best1(N)` is the minimum of the
native RMSD over the first N of the score order: it says the candidate IS there, not that anything
native-free can find it.  No deployable parameter is chosen anywhere in this file.

    python s31/s31_C_widen.py run
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

from s12 import instrument as I                    # noqa: E402
from s24 import stats_lib as ST                    # noqa: E402
from s31 import s31_C_cache as CA                  # noqa: E402

RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "s31_C_widen.json")
N_NULL = 20000
NARROW, WIDE = 128, 500          # the deployed register, and the whole shipped pool


def chain_rows():
    """s29 lane O's already-built chains, keyed item -> pdb -> chain RMSD.  Nothing is rebuilt."""
    import glob
    out = {}
    for f in glob.glob(os.path.join(ROOT, "s29", "results", "s29_O_chain_rows*.jsonl")):
        for line in open(f):
            r = json.loads(line)
            out.setdefault(r["item"], {})[r["pdb"]] = float(r["rmsd_chain"])
    return out


#: SIGN CORRECTION, made before any number left this lane.  The first version of this file defined
#: `g = best1(128) - best1(500)`, which is >= 0 by construction (a running minimum cannot rise), and
#: then read the falsifiers and the null's 2.5th percentile as though more negative were better --
#: the WRONG TAIL of the right null, which is S31 contract rule 8's exact failure mode.  The
#: registered quantity is the EFFECT OF WIDENING, `eff = best1(500) - best1(128) <= 0`, matching the
#: sign of S30-L11's published -1.9004.  Everything below is in `eff`.
def stratum_null(g, idx_size, n_draw, seed_parts):
    """Distribution of the mean of `g` over a random `idx_size`-subset.  The matched null for a
    stratum mean (contract rule 8: the statistic and its null are the same object)."""
    rng = ST._rng(*seed_parts)
    n = len(g)
    draws = np.array([g[rng.choice(n, idx_size, replace=False)].mean() for _ in range(n_draw)])
    return draws


def main():
    pdbs = CA.all_pdbs()
    rows = []
    for p in pdbs:
        z = CA.load(p)
        order = np.asarray(z["order"], int)
        rr = np.asarray(z["rr"], float)                                  # ORACLE
        run = np.minimum.accumulate(rr[order])                           # best1(N) for N = 1..k
        best_rank = int(np.argmin(rr[order])) + 1                        # 1-based rank of the argmin
        rows.append(dict(
            pdb=p, n=int(z["n"]), fold=int(z["fold"]),
            best1_75=float(run[74]), best1_128=float(run[NARROW - 1]),
            best1_256=float(run[255]), best1_500=float(run[WIDE - 1]),
            pool_mean=float(rr.mean()), best_rank=best_rank,
            curve={str(N): float(run[N - 1]) for N in (1, 2, 4, 8, 16, 32, 64, 128, 256, 500)},
        ))

    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    b128 = np.array([r["best1_128"] for r in rows])
    b500 = np.array([r["best1_500"] for r in rows])
    b75 = np.array([r["best1_75"] for r in rows])
    pmean = np.array([r["pool_mean"] for r in rows])
    rank = np.array([r["best_rank"] for r in rows])
    g = b500 - b128                                # EFFECT OF WIDENING, <= 0; ORACLE

    # ------------------------------------------------------------------ strata
    fail18 = np.array([p in I.FAIL18 for p in pdbs])
    worst_pm = np.zeros(len(rows), bool); worst_pm[np.argsort(-pmean)[:18]] = True
    worst_bp = np.zeros(len(rows), bool); worst_bp[np.argsort(-b500)[:18]] = True
    defn18 = np.zeros(len(rows), bool); defn18[np.argsort(-(b75 - b500))[:18]] = True
    strata = {"FAIL18": fail18, "worst18_poolmean": worst_pm,
              "worst18_bestpool": worst_bp, "defn18": defn18}

    null = stratum_null(g, 18, N_NULL, ("s31C", "widen", "null18"))
    out = {"n": len(rows), "basis": "CA point cloud unless a key says chain",
           "ORACLE": True, "deployable": False,
           "all126": {"g": float(g.mean()), "se": float(g.std(ddof=1) / np.sqrt(len(g))),
                      "best1_128": float(b128.mean()), "best1_500": float(b500.mean())},
           "null18": {"mean": float(null.mean()), "sd": float(null.std(ddof=1)),
                      "p2.5": float(np.percentile(null, 2.5)),
                      "p97.5": float(np.percentile(null, 97.5)), "n_draw": N_NULL},
           "strata": {}, "overlap": {}}

    for name, m in strata.items():
        gm = float(g[m].mean())
        # one-sided p: how often does a random 18 reach a mean at least this negative
        p = float((null <= gm).mean())      # lower tail: more negative = widening helps
        out["strata"][name] = {
            "n": int(m.sum()), "g": gm,
            "g_se": float(g[m].std(ddof=1) / np.sqrt(m.sum())),
            "best1_128": float(b128[m].mean()), "best1_500": float(b500[m].mean()),
            "g_norm": float((g[m] / b128[m]).mean()),
            "p_vs_null18": p,
            "clears_null_2.5pct": bool(gm < out["null18"]["p2.5"]),
            "frac_rank_gt_128": float((rank[m] > NARROW).mean()),
            "median_rank": float(np.median(rank[m])),
            "members": [pdbs[i] for i in np.where(m)[0]],
        }
    out["strata"]["other108"] = {
        "n": int((~fail18).sum()), "g": float(g[~fail18].mean()),
        "g_se": float(g[~fail18].std(ddof=1) / np.sqrt((~fail18).sum())),
        "best1_128": float(b128[~fail18].mean()), "best1_500": float(b500[~fail18].mean()),
        "frac_rank_gt_128": float((rank[~fail18] > NARROW).mean()),
        "median_rank": float(np.median(rank[~fail18])),
    }
    out["all126_rank"] = {"frac_rank_gt_128": float((rank > NARROW).mean()),
                          "median_rank": float(np.median(rank)),
                          "uninformative_null_frac": float((WIDE - NARROW) / WIDE)}

    for a in strata:
        for b in strata:
            if a < b:
                out["overlap"]["%s&%s" % (a, b)] = int((strata[a] & strata[b]).sum())

    # ------------------------------------------------------------------ S30's published contrast
    d = g[fail18].mean() - g[~fail18].mean()
    se = float(np.sqrt(g[fail18].var(ddof=1) / fail18.sum() + g[~fail18].var(ddof=1) / (~fail18).sum()))
    out["S30_published_contrast"] = {"fail18_minus_108": float(d), "se": se,
                                     "x_mde": float(abs(d) / (2.8016 * se))}

    # ------------------------------------------------------------------ level control
    A = np.vstack([np.ones(len(g)), b128]).T
    beta, *_ = np.linalg.lstsq(A, g, rcond=None)
    resid = g - A @ beta
    out["level_control"] = {
        "slope_on_best1_128": float(beta[1]), "intercept": float(beta[0]),
        "r2": float(1 - resid.var() / g.var()),
        "residual_by_stratum": {k: float(resid[m].mean()) for k, m in strata.items()},
        "residual_other108": float(resid[~fail18].mean()),
        "note": "g is mechanically bounded by best1_128; the residual is the part not explained "
                "by how bad the 128-window already is",
    }

    # ------------------------------------------------------------------ fold-clustered CI on the
    # headline stratum, computed on the paired arms rather than on g alone
    m = worst_pm
    out["headline_fold"] = ST.compare(b500[m], b128[m], folds=folds[m],
                                      names=[pdbs[i] for i in np.where(m)[0]],
                                      label="widen_worst18_poolmean",
                                      seed_parts=("s31C", "widen"))
    m = fail18
    out["fail18_fold"] = ST.compare(b500[m], b128[m], folds=folds[m],
                                    names=[pdbs[i] for i in np.where(m)[0]],
                                    label="widen_fail18", seed_parts=("s31C", "widen"))

    # ------------------------------------------------------------------ the BUILT CHAIN arm
    ch = chain_rows()
    if "best1_top128" in ch and "best1_pool" in ch:
        c128 = np.array([ch["best1_top128"][p] for p in pdbs])
        c500 = np.array([ch["best1_pool"][p] for p in pdbs])
        gc = c500 - c128                               # EFFECT OF WIDENING on the chain
        nullc = stratum_null(gc, 18, N_NULL, ("s31C", "widen", "nullchain"))
        out["chain"] = {"basis": "BUILT CHAIN (s29 lane O rows; nothing rebuilt)",
                        "all126": float(gc.mean()),
                        "null18_p2.5": float(np.percentile(nullc, 2.5))}
        for name, mm in list(strata.items()) + [("other108", ~fail18)]:
            out["chain"][name] = {"g": float(gc[mm].mean()),
                                  "p_vs_null18": float((nullc <= gc[mm].mean()).mean())}

    # ------------------------------------------------------------------ verdict
    h = out["strata"]["worst18_poolmean"]
    if h["g"] <= -1.00 and h["clears_null_2.5pct"]:
        v = "F-C3a FIRES -- REVIVE"
    elif h["g"] >= -0.50 or not h["clears_null_2.5pct"]:
        v = "F-C3b FIRES -- CLOSE"
    else:
        v = "PARTIAL -- survives at reduced size"
    out["verdict"] = v
    out["definitional_diagnostic"] = {
        "defn18_g": out["strata"]["defn18"]["g"], "fail18_g": out["strata"]["FAIL18"]["g"],
        "defn18_at_least_as_large": bool(abs(out["strata"]["defn18"]["g"]) >=
                                         abs(out["strata"]["FAIL18"]["g"])),
    }

    os.makedirs(RESULTS, exist_ok=True)
    out["provenance"] = ST.provenance(__file__)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    with open(os.path.join(RESULTS, "s31_C_widen_rows.jsonl"), "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

    # ------------------------------------------------------------------ render
    print("C3 -- the 128->512 widening.  EVERY ROW ORACLE / NOT DEPLOYABLE.  CA point cloud.\n")
    print("  g = best1(128) - best1(500);  more negative = widening helps more\n")
    print("  %-20s %4s %9s %9s %9s %9s %8s %8s" %
          ("stratum", "n", "g", "best1_128", "best1_500", "g/b128", "p_null", "rk>128"))
    for k in ("FAIL18", "other108", "worst18_poolmean", "worst18_bestpool", "defn18"):
        s = out["strata"][k]
        print("  %-20s %4d %9.4f %9.4f %9.4f %9.4f %8.4f %8.3f" %
              (k, s["n"], s["g"], s["best1_128"], s["best1_500"], s.get("g_norm", float("nan")),
               s.get("p_vs_null18", float("nan")), s["frac_rank_gt_128"]))
    print("  %-20s %4d %9.4f %9.4f %9.4f" %
          ("all126", out["n"], out["all126"]["g"], out["all126"]["best1_128"],
           out["all126"]["best1_500"]))
    print("\n  random-18 null on g: mean %.4f sd %.4f  2.5pct %.4f" %
          (out["null18"]["mean"], out["null18"]["sd"], out["null18"]["p2.5"]))
    print("  S30 published contrast FAIL18-108: %.4f  SE %.4f  %.2fx MDE" %
          (out["S30_published_contrast"]["fail18_minus_108"],
           out["S30_published_contrast"]["se"], out["S30_published_contrast"]["x_mde"]))
    print("  level control: g ~ best1_128 slope %.4f R2 %.4f; residual by stratum %s" %
          (out["level_control"]["slope_on_best1_128"], out["level_control"]["r2"],
           {k: round(v, 4) for k, v in out["level_control"]["residual_by_stratum"].items()}))
    print("  stratum overlaps: %s" % out["overlap"])
    if "chain" in out:
        print("\n  BUILT CHAIN (s29 rows): all126 %.4f  FAIL18 %.4f  worst18_poolmean %.4f  "
              "worst18_bestpool %.4f  other108 %.4f" %
              (out["chain"]["all126"], out["chain"]["FAIL18"]["g"],
               out["chain"]["worst18_poolmean"]["g"], out["chain"]["worst18_bestpool"]["g"],
               out["chain"]["other108"]["g"]))
    print("\n  VERDICT: %s" % v)
    print("  wrote %s" % OUT)


if __name__ == "__main__":
    main()
