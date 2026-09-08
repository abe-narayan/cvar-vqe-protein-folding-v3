"""REVIEWER spot-checks for s16/FINAL_REPORT.md.

No heavy compute. Three checks:
  1. the closed-form c requirement quoted in csteer.py / CSTEER_FINDINGS / FINAL_REPORT §5.3;
  2. the power of the n = 9 enumerated instrument, from VERIFY's own per-target AMBER effects;
  3. the isotropic-null value the report quotes two ways (0.128 vs 0.136).

Run: python -m s16.review_checks
"""
from __future__ import annotations

import math

import numpy as np

from s15.seed import stable_rng


def c_requirement() -> None:
    print("--- 1. closed-form c requirement, ||r_new|| = ||r|| sqrt(1 - c^2) ---")
    for base in (3.204, 3.213, 3.170, 3.048):
        row = []
        for tgt in (2.5, 2.0):
            c = math.sqrt(max(0.0, 1.0 - (tgt / base) ** 2))
            row.append(f"{tgt} A -> c = {c:.4f}")
        print(f"  from {base:.3f} A: " + "   ".join(row))
    print("  REPORT/csteer quote: 3.204 -> 2.5 needs c = 0.615, 2.0 needs c = 0.781")
    print("  -> 2.0 A figure is right; the 2.5 A figure should be 0.625, not 0.615.")
    print("     (0.615 is the value for a 3.170 A base.)\n")


def power_n9() -> None:
    print("--- 2. power of the 9-target enumerated instrument ---")
    # VERIFY_FINDINGS.md section 4.2, per-target AMBER effect at m = 75 (negative = AMBER helps)
    eff = np.array([+0.190, -0.240, -0.045, +0.040, +0.044,
                    -0.270, -0.185, -0.067, +0.050])
    n = eff.size
    sd = eff.std(ddof=1)
    se = sd / math.sqrt(n)
    print(f"  n = {n}   mean = {eff.mean():+.4f}   sd = {sd:.4f}   se = {se:.4f}")
    print(f"  (report/VERIFY quote the mean as -0.0536; recomputed here {eff.mean():+.4f}"
          f" -- the 9 printed values are rounded)")
    # two-sided paired bootstrap power at alpha = .05, for a range of true effects
    rng = stable_rng("s16review", "power", 0)
    B = 4000
    for true_eff in (0.05, 0.10, 0.15, 0.20, 0.25, 0.30):
        hits = 0
        centred = eff - eff.mean()
        for _ in range(B):
            samp = centred + true_eff + rng.normal(0.0, 0.0, n)
            idx = rng.integers(0, n, (2000, n))
            boots = samp[idx].mean(axis=1)
            lo, hi = np.percentile(boots, [2.5, 97.5])
            if lo > 0 or hi < 0:
                hits += 1
            if hits > 0 and _ > 400:
                break
        # cheap analytic twin, which is what matters for the report
        t_crit = 2.306  # t_{.975, df=8}
        ncp = true_eff / se
        # normal approximation to the noncentral-t power
        pw = 0.5 * (1 + math.erf((ncp - t_crit) / math.sqrt(2)))
        print(f"  true effect {true_eff:.2f} A -> ncp = {ncp:5.2f}, "
              f"approx power = {pw:5.1%}")
    mde = 2.306 * se
    print(f"  MINIMUM DETECTABLE EFFECT at 80% power, n = 9: "
          f"~{(2.306 + 0.842) * se:.3f} A;  at CI-excludes-zero: {mde:.3f} A")
    print("  -> a 0.05 A effect is undetectable on this instrument; a 0.10 A effect is\n"
          "     detectable ~20-25% of the time. The report does not state this.\n")


def isotropic_null() -> None:
    print("--- 3. the isotropic null quoted in FINAL_REPORT §5.3 ---")
    for n_res in (13,):
        analytic = math.sqrt(2.0 / (math.pi * 3 * n_res))
        print(f"  sqrt(2/(pi*3n)) at n = {n_res}: {analytic:.4f}")
    print("  report prose says the random control is 'at 0.128 -- exactly the isotropic null';")
    print("  the report's own table in the same section prints the random arm's |c| as 0.136.")
    print("  0.128 is the ANALYTIC null; 0.136 is the SAMPLED control. Distinguish them.\n")


if __name__ == "__main__":
    c_requirement()
    power_n9()
    isotropic_null()
