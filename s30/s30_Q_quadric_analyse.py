#!/usr/bin/env python
"""s30/s30_Q_quadric_analyse.py -- lane Q: does the QUADRIC (second-moment) reachable-tail class
beat the HALFSPACE class? Reads s30/results/s30_Q_quadric_rows.jsonl. POINT CLOUD, ALL ORACLE.

The readout is held fixed (uniform average of the selected prefix) so the contrast isolates WHICH
SET -- the only thing the cut class changes. Order-statistic discipline is applied to the sampled
direction search, which is where the apparent ceiling comes from.
"""
from __future__ import annotations
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
if ROOT not in sys.path: sys.path.insert(0, ROOT)
from s12 import instrument as I      # noqa: E402
from s24 import stats_lib as ST      # noqa: E402

ROWS = os.path.join(HERE, "results", "s30_Q_quadric_rows.jsonl")
OUT = os.path.join(HERE, "results", "s30_Q_quadric.json")
K_CURVE = (1, 2, 4, 8, 16, 32, 64, 128)


def main():
    rows = {}
    for L in open(ROWS):
        L = L.strip()
        if L:
            r = json.loads(L); rows[r["pdb"]] = r
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in rows]
    n = len(pdbs); fm = np.array([p in I.FAIL18 for p in pdbs])
    g = lambda f: np.array([f(rows[p]) for p in pdbs], float)
    out = dict(n=n, provenance=ST.provenance(__file__))
    exo = g(lambda r: r["classes"]["exogenous"]["best"])
    exo75 = g(lambda r: r["classes"]["exogenous"]["at_m75"])
    osort = g(lambda r: r["classes"]["oracle_sort"]["best"])
    print("== QUADRIC vs HALFSPACE reachable-tail class, n=%d. POINT CLOUD. ALL ORACLE." % n)
    print("   deployed exogenous prefix, ORACLE m  %.4f   (shipped m=75: %.4f)" % (exo.mean(), exo75.mean()))
    print("   ORACLE-sorted prefix (upper ref)     %.4f\n" % osort.mean())
    out.update(exogenous_oracle_m=float(exo.mean()), exogenous_m75=float(exo75.mean()),
               oracle_sort=float(osort.mean()))
    print("   best-of-K over sampled directions (ORDER STATISTIC -- read the curve):")
    print("   %-6s %10s %10s %10s %10s" % ("K", "LIN best", "QUAD best", "LIN @m75", "QUAD @m75"))
    out["curveK"] = {}
    for i, K in enumerate(K_CURVE):
        v = [g(lambda r: r["classes"][c]["%s" % k][i]).mean()
             for c, k in (("linear", "curve_bestK"), ("quadric", "curve_bestK"),
                          ("linear", "at_m75_curveK"), ("quadric", "at_m75_curveK"))]
        out["curveK"][K] = dict(lin_best=v[0], quad_best=v[1], lin_m75=v[2], quad_m75=v[3])
        print("   %-6d %10.4f %10.4f %10.4f %10.4f" % (K, *v))
    LB = g(lambda r: r["classes"]["linear"]["best"]); QB = g(lambda r: r["classes"]["quadric"]["best"])
    L7 = g(lambda r: r["classes"]["linear"]["at_m75"]); Q7 = g(lambda r: r["classes"]["quadric"]["at_m75"])
    for lbl, a, b in (("ORACLE m", QB, LB), ("fixed m=75", Q7, L7)):
        d = a - b; se = d.std(ddof=1) / math.sqrt(n)
        out["quad_minus_lin_%s" % lbl.split()[0]] = dict(mean=float(d.mean()), se=float(se),
                                                         mde=float(abs(d.mean()) / (2.8016 * se)),
                                                         wins=int((d < 0).sum()), losses=int((d > 0).sum()))
        print("\n   QUADRIC minus LINEAR at %-11s %+.4f  SE %.4f  -> %.2fx MDE  (W/L %d/%d)"
              % (lbl, d.mean(), se, abs(d.mean()) / (2.8016 * se), (d < 0).sum(), (d > 0).sum()))
    print("\n   best_of_k_within on the sampled direction search (128 directions):")
    out["order_stat"] = {}
    for c in ("linear", "quadric"):
        M = np.array([rows[p]["classes"][c]["per_dir_best"] for p in pdbs])
        bk = ST.best_of_k_within(M, seed_parts=("s30Q", c)); out["order_stat"][c] = bk
        print("     %-8s observed %+.4f  across-target null %+.4f (%.0f%% accounted)  "
              "split-half %+.4f (%.0f%%)  k_eff %.1f"
              % (c, bk["observed_gain"], bk["null_across_targets"], 100 * bk["share_accounted"],
                 bk["split_half"], 100 * bk["split_half_frac"], bk["k_eff"]))
    ml = np.concatenate([rows[p]["classes"]["linear"]["per_dir_m"] for p in pdbs])
    out["lin_best_m_median"] = float(np.median(ml)); out["lin_frac_m_ge50"] = float((ml >= 50).mean())
    print("\n   best prefix for a LINEAR direction: median m %d, %.0f%% at m<=2, %.0f%% at m>=50"
          % (int(np.median(ml)), 100 * (ml <= 2).mean(), 100 * (ml >= 50).mean()))
    print("\n   NATIVE-FREE RULE directions (the only DEPLOYABLE rows):")
    out["rules"] = {}
    res = []
    for nm in sorted(rows[pdbs[0]]["rules"]):
        b = g(lambda r: r["rules"][nm]["best"]); a75 = g(lambda r: r["rules"][nm]["at_m75"])
        out["rules"][nm] = dict(oracle_m=float(b.mean()), m75=float(a75.mean()))
        res.append((a75.mean(), b.mean(), nm))
    for a75, b, nm in sorted(res)[:8]:
        print("     %-18s ORACLE-m %.4f   at m=75 %.4f" % (nm, b, a75))
    print("\n   FAIL18 / 108 at fixed m=75:")
    out["strata_m75"] = {}
    for nm, v in (("exogenous", exo75), ("LIN best-of-128", L7), ("QUAD best-of-128", Q7)):
        out["strata_m75"][nm] = dict(all=float(v.mean()), fail18=float(v[fm].mean()),
                                     other108=float(v[~fm].mean()))
        print("     %-18s all %.4f  FAIL18 %.4f  other108 %.4f" % (nm, v.mean(), v[fm].mean(), v[~fm].mean()))
    ST.save_atomic(OUT, out, module_file=__file__)
    print("\nwrote", OUT)
    return out


if __name__ == "__main__":
    main()
