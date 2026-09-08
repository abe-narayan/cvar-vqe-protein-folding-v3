"""AGENT D, Sprint 19 -- T3 report: how much of the tangent story is a model and how much a reference.

The coordinator asked for exactly two things before promoting the design consequence past SUPPORTED:

  1. does the per-target in-tangent EXCESS over the analytic null predict the per-target
     Sprint-18 gap, or does the statistic only order the ARMS?
  2. an honest statement of the linearisation error -- the fit moves several radians in torsion
     space, which is not a neighbourhood, so how well does a tangent projection at the START
     predict what the fit actually REALISES (AGENT A's kappa)?

Both are answered here against `s19/results/D_T3_tangent/tangent.json` and `s19/results/a_coh.json`.

Run:  python -m s19.d_tangent_report
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
from s15 import seed as SD                  # noqa: E402

ARMS = ["real", "signflip", "shuffled", "iso", "coherent"]


def main():
    T = json.load(open(os.path.join(HERE, "results", "D_T3_tangent", "tangent.json")))
    A = json.load(open(os.path.join(HERE, "results", "a_coh.json")))
    assert T.get("complete") and A.get("complete"), "refusing to report a partial file"
    a = {r["pdb"]: r for r in A["rows"]}
    rows = [r for r in T["rows"] if r["pdb"] in a]
    rng = SD.stable_rng("D", "T3rep")
    from scipy.stats import spearmanr

    def boot(v, B=4000):
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        b = v[rng.integers(0, len(v), size=(B, len(v)))].mean(1)
        return float(v.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))

    g = lambda k: np.asarray([r[k] for r in rows], float)          # noqa: E731
    ga = lambda k: np.asarray([a[r["pdb"]][k] for r in rows], float)   # noqa: E731

    print(f"\nn = {len(rows)} targets\n")
    print("== 2. LINEARISATION ERROR: the tangent fraction at the START vs what the fit REALISES")
    print(f"  {'arm':<12}{'in-tangent @start':>19}{'A kappa (realised)':>21}"
          f"{'diff':>9}{'per-target rho':>16}")
    for arm in ARMS:
        f0, kp = g("fin_" + arm), ga(arm + "_kappa")
        print(f"  {arm:<12}{f0.mean():>19.3f}{kp.mean():>21.3f}{kp.mean()-f0.mean():>+9.3f}"
              f"{spearmanr(f0, kp).statistic:>+16.3f}")
    f0 = np.concatenate([g("fin_" + x) for x in ARMS])
    kp = np.concatenate([ga(x + "_kappa") for x in ARMS])
    print(f"  POOLED across arms x targets: rho = {spearmanr(f0, kp).statistic:+.3f}, "
          f"pearson {np.corrcoef(f0, kp)[0, 1]:+.3f}, mean|diff| {np.abs(f0-kp).mean():.3f}")
    am = np.array([g("fin_" + x).mean() for x in ARMS])
    ak = np.array([ga(x + "_kappa").mean() for x in ARMS])
    print(f"  across the {len(ARMS)} ARM MEANS ONLY:   pearson {np.corrcoef(am, ak)[0, 1]:+.3f}")

    print("\n== the analytic null, and whether the incoherent arms sit on it")
    m, lo, hi = boot(g("null_DP"))
    print(f"  D/P = (2n-2)/npairs                {m:.3f}  [{lo:.3f},{hi:.3f}]")
    m, lo, hi = boot(g("rank_Jw") / g("npairs"))
    print(f"  rank(J_w)/npairs (the TRUE ratio)  {m:.3f}  [{lo:.3f},{hi:.3f}]")
    for arm in ARMS + ["random"]:
        m, lo, hi = boot(g("fin_" + arm))
        print(f"  in-tangent, {arm:<22}{m:.3f}  [{lo:.3f},{hi:.3f}]")
    m, lo, hi = boot(g("fin_real") - g("fin_random"))
    print(f"  real - empirical random direction  {m:+.3f} [{lo:+.3f},{hi:+.3f}]  "
          f"W/L {int((g('fin_real')>g('fin_random')).sum())}/"
          f"{int((g('fin_real')<g('fin_random')).sum())}")

    print("\n== 1. DOES THE EXCESS PREDICT THE PER-TARGET GAP, or does it only order the ARMS?")
    gap = ga("real") - ga("signflip")             # the Sprint-18 gap, per target
    exc = g("fin_real") - g("fin_random")
    ex2 = g("fin_real") - g("null_DP")
    print(f"  gap (real - signflip) mean {gap.mean():+.3f}, sd {gap.std():.3f}")
    for lab, v in [("in-tangent excess over empirical random", exc),
                   ("in-tangent excess over D/P", ex2),
                   ("in-tangent fraction, raw", g("fin_real")),
                   ("A's kappa(real), raw", ga("real_kappa")),
                   ("A's kappa(real) - kappa(signflip)", ga("real_kappa") - ga("signflip_kappa"))]:
        r_ = spearmanr(v, gap)
        print(f"  rho({lab:<40}, gap) = {r_.statistic:+.3f}  p = {r_.pvalue:.3g}")
    # a bootstrap CI on the headline correlation
    n = len(gap)
    bs = np.array([spearmanr(exc[k], gap[k]).statistic
                   for k in rng.integers(0, n, size=(2000, n))])
    print(f"  95% CI on rho(excess over random, gap) = "
          f"[{np.percentile(bs,2.5):+.3f}, {np.percentile(bs,97.5):+.3f}]")

    json.dump({"n": len(rows)}, open(os.path.join(HERE, "results", "D_T3_tangent",
                                                  "report_meta.json"), "w"))


if __name__ == "__main__":
    main()
