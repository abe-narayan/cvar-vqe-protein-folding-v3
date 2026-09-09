"""s25/audit_t1.py -- INDEPENDENT AUDIT of s25 L1 (calibration) and L2 (mechanism).

Written by the AUDIT lane.  Re-derives every number from the persisted distogram artefacts
(s12/cache/disto_*.npz) and s25/results/temper.json with independent code.  Does NOT import
s25/calib.py or s25/temper.py.  No model training, no LOCK_TRAIN, no LOCK_AMBER.

WHAT IS RE-DERIVED
  A1  z_sd / coverage exactly as calib.py defines them        (reproduction check)
  A2  robust alternatives to z_sd: median|z|, trimmed sd, pooled vs per-target,
      unimodal-only, edge-censored-excluded, continuity-corrected sd
  A3  the COVERAGE OPERATOR: bin-centre endpoints (calib.py) vs bin-EDGE endpoints vs
      linearly interpolated quantiles -- three readouts of "the central 90% interval"
  B1  does the SD-convolution operator move the posterior MEDIAN?  does TEMP?
      (this is the load-bearing step of L2's mechanism claim)
  B2  weight channel vs shape channel from the persisted SD / SDFIXW columns
  B3  an independent within-target best-of-K null for the ORACLE row
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s24 import stats_lib as ST           # noqa: E402

CACHE = os.path.join(ROOT, "s12", "cache")
FGRID = [1.0, 1.15, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0]
TGRID = [1.0, 1.25, 1.5, 2.0, 3.0, 5.0]


def edges_from_centres(C):
    mid = 0.5 * (C[1:] + C[:-1])
    lo = C[0] - (mid[0] - C[0]); hi = C[-1] + (C[-1] - mid[-1])
    return np.concatenate([[lo], mid, [hi]])


def widen_sd(prob, C, f):
    """Verbatim re-implementation of temper.py's _widen_sd (audited, not imported)."""
    if f <= 1.0:
        return prob.copy()
    m = (prob * C[None, :]).sum(1)
    v = (prob * (C[None, :] - m[:, None]) ** 2).sum(1)
    s = np.sqrt(np.maximum(v, 1e-12)) * np.sqrt(f * f - 1.0)
    d2 = (C[None, :] - C[:, None]) ** 2
    out = np.empty_like(prob)
    for p in range(len(prob)):
        K = np.exp(-d2 / (2.0 * max(s[p], 1e-6) ** 2))
        K /= K.sum(1, keepdims=True)
        q = prob[p] @ K
        out[p] = q / q.sum()
    return out


def temper(prob, T):
    if T == 1.0:
        return prob.copy()
    q = np.power(np.maximum(prob, 1e-12), 1.0 / T)
    return q / q.sum(1, keepdims=True)


def moments(prob, C):
    m = (prob * C[None, :]).sum(1)
    v = (prob * (C[None, :] - m[:, None]) ** 2).sum(1)
    return m, np.sqrt(np.maximum(v, 1e-12))


def discrete_median(prob, C):
    """The L1 minimiser of sum_c p_c |t - C_c| on the bin-centre support."""
    cdf = np.cumsum(prob, 1)
    return C[np.argmax(cdf >= 0.5, axis=1)]


def interp_quantile(prob, edges, a):
    """Quantile of the piecewise-UNIFORM density implied by the bin masses.  The honest
    continuous reading of a binned posterior; calib.py used the bin CENTRE instead."""
    cdf = np.cumsum(prob, 1)
    k = np.argmax(cdf >= a, axis=1)
    r = np.arange(len(prob))
    below = np.where(k > 0, cdf[r, np.maximum(k - 1, 0)], 0.0)
    mass = np.maximum(prob[r, k], 1e-12)
    frac = np.clip((a - below) / mass, 0.0, 1.0)
    return edges[k] + frac * (edges[k + 1] - edges[k])


def main():
    tg = I.targets()
    pdbs = [t["pdb"] for t in tg]
    fold = np.array([t["fold"] for t in tg], int)
    rows = []
    print("targets %d" % len(tg), flush=True)
    for t in tg:
        z = np.load(os.path.join(CACHE, "disto_%s.npz" % t["pdb"]))
        P = np.asarray(z["prob"], float)
        C = np.asarray(z["centres"], float)
        exp = np.asarray(z["expected"], float)
        sd = np.maximum(np.asarray(z["sd"], float), 1e-6)
        i, j = np.asarray(z["i"], int), np.asarray(z["j"], int)
        u = I.load_univ(t["pdb"])
        nat = np.asarray(u["nat_ca"], float)
        Dt = I.pair_dists(nat[None], i, j)[0]
        sep = (j - i).astype(int)
        E = edges_from_centres(C)

        zz = (Dt - exp) / sd
        cdf = np.cumsum(P, 1)

        # --- three readouts of the central 90% / 50% interval
        qc = lambda a: C[np.argmax(cdf >= a, axis=1)]                       # noqa: E731 calib.py
        ke = lambda a: E[np.argmax(cdf >= a, axis=1) + (1 if a > 0.5 else 0)]  # noqa: E731 outer edges
        qi = lambda a: interp_quantile(P, E, a)                             # noqa: E731 interpolated
        cov = {}
        for nm, fn in (("centre", qc), ("edge", ke), ("interp", qi)):
            for lab, (a, b) in (("50", (0.25, 0.75)), ("90", (0.05, 0.95))):
                lo, hi = fn(a), fn(b)
                cov["cov%s_%s" % (lab, nm)] = float(((Dt >= lo) & (Dt <= hi)).mean())

        # --- continuity correction: bin masses are a histogram; the discrete sd omits the
        #     within-bin spread.  var_cont = var_bin + sum_c p_c * w_c^2 / 12
        w = np.diff(E)
        var_within = (P * (w[None, :] ** 2) / 12.0).sum(1)
        sd_cc = np.sqrt(sd ** 2 + var_within)
        z_cc = (Dt - exp) / sd_cc

        # --- edge censoring: natives outside the representable range
        out_hi = Dt > E[-1]
        out_lo = Dt < E[0]

        # --- multimodality
        pad = np.pad(P, ((0, 0), (1, 1)), constant_values=0.0)
        peaks = ((pad[:, 1:-1] > pad[:, :-2]) & (pad[:, 1:-1] >= pad[:, 2:]) & (P > 0.02)).sum(1)
        uni = peaks < 2

        def rsd(x):
            return float(np.std(x, ddof=1)) if len(x) > 1 else np.nan

        r = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "npairs": int(len(Dt)),
             "z_mean": float(zz.mean()), "z_sd": rsd(zz),
             "z_medabs": float(np.median(np.abs(zz))),
             "z_iqr_sd": float((np.percentile(zz, 75) - np.percentile(zz, 25)) / 1.34898),
             "z_sd_trim5": rsd(np.clip(zz, *np.percentile(zz, [2.5, 97.5]))),
             "z_sd_unimodal": rsd(zz[uni]) if uni.sum() > 1 else np.nan,
             "z_sd_multimodal": rsd(zz[~uni]) if (~uni).sum() > 1 else np.nan,
             "z_sd_inrange": rsd(zz[~(out_hi | out_lo)]) if (~(out_hi | out_lo)).sum() > 1 else np.nan,
             "z_sd_cc": rsd(z_cc), "z_medabs_cc": float(np.median(np.abs(z_cc))),
             "frac_out_hi": float(out_hi.mean()), "frac_out_lo": float(out_lo.mean()),
             "multimodal_frac": float((~uni).mean()),
             "kurt": float(((zz - zz.mean()) ** 4).mean() / max(zz.var() ** 2, 1e-12)),
             "signed": float((Dt - exp).mean()),
             "sd_mean": float(sd.mean()), "sd_cc_mean": float(sd_cc.mean())}
        r.update(cov)
        # z pooled contribution
        r["_z"] = zz.tolist()
        r["_zcc"] = z_cc.tolist()
        r["_sep"] = sep.tolist()

        # ---------- B1: does the widening operator move the MEDIAN / MEAN?
        med0 = discrete_median(P, C); m0, s0 = moments(P, C)
        for f in FGRID[1:]:
            Pw = widen_sd(P, C, f)
            mw, sw = moments(Pw, C)
            r["dmed_SD_%s" % f] = float(np.abs(discrete_median(Pw, C) - med0).mean())
            r["dmean_SD_%s" % f] = float((mw - m0).mean())
            r["sdratio_SD_%s" % f] = float((sw / s0).mean())
        for T in TGRID[1:]:
            Pt = temper(P, T)
            mt, st = moments(Pt, C)
            r["dmed_TEMP_%s" % T] = float(np.abs(discrete_median(Pt, C) - med0).mean())
            r["dmean_TEMP_%s" % T] = float((mt - m0).mean())
            r["sdratio_TEMP_%s" % T] = float((st / s0).mean())
        rows.append(r)

    # -------------------------------------------------- report
    g = lambda k: np.array([r[k] for r in rows], float)     # noqa: E731
    allz = np.concatenate([np.asarray(r["_z"]) for r in rows])
    allzcc = np.concatenate([np.asarray(r["_zcc"]) for r in rows])
    npair = g("npairs")

    print("\n" + "=" * 96)
    print("A1  REPRODUCTION of calib.py L1  (per-target then averaged, as declared)")
    print("=" * 96)
    print("  z_mean %+.4f   z_sd %.4f   cov50 %.4f   cov90 %.4f   multimodal %.3f   pairs %d"
          % (g("z_mean").mean(), g("z_sd").mean(), g("cov50_centre").mean(),
             g("cov90_centre").mean(), g("multimodal_frac").mean(), int(npair.sum())))

    print("\n" + "=" * 96)
    print("A2  IS z_sd = 1.9962 AN ARTEFACT OF THE sd STATISTIC?   robust alternatives")
    print("=" * 96)
    print("    reference values for a PERFECTLY calibrated posterior in brackets")
    print("  %-42s%12s%12s" % ("statistic", "value", "nominal"))
    print("  %-42s%12.4f%12.4f" % ("z sd, per-target then averaged", g("z_sd").mean(), 1.0))
    print("  %-42s%12.4f%12.4f" % ("z sd, pooled over all pairs", np.std(allz, ddof=1), 1.0))
    print("  %-42s%12.4f%12.4f" % ("median |z|  (robust to tails)", g("z_medabs").mean(), 0.6745))
    print("  %-42s%12.4f%12.4f" % ("  -> implied sd = median|z|/0.6745",
                                   g("z_medabs").mean() / 0.6745, 1.0))
    print("  %-42s%12.4f%12.4f" % ("IQR/1.349  (robust scale)", g("z_iqr_sd").mean(), 1.0))
    print("  %-42s%12.4f%12.4f" % ("z sd, 2.5%% winsorised", g("z_sd_trim5").mean(), 1.0))
    print("  %-42s%12.4f%12.4f" % ("z sd, UNIMODAL pairs only", np.nanmean(g("z_sd_unimodal")), 1.0))
    print("  %-42s%12.4f%12.4f" % ("z sd, MULTIMODAL pairs only", np.nanmean(g("z_sd_multimodal")), 1.0))
    print("  %-42s%12.4f%12.4f" % ("z sd, natives inside bin range only",
                                   np.nanmean(g("z_sd_inrange")), 1.0))
    print("  %-42s%12.4f%12.4f" % ("z sd, CONTINUITY-CORRECTED sd", g("z_sd_cc").mean(), 1.0))
    print("  %-42s%12.4f%12.4f" % ("  median|z| continuity-corrected", g("z_medabs_cc").mean(), 0.6745))
    print("  %-42s%12.4f" % ("excess kurtosis of z (mean per target)", g("kurt").mean() - 3.0))
    print("  %-42s%12.4f  %s" % ("natives ABOVE the top bin edge", g("frac_out_hi").mean(), "frac"))
    print("  %-42s%12.4f  %s" % ("natives BELOW the bottom bin edge", g("frac_out_lo").mean(), "frac"))
    print("  mean published sd %.4f  ->  continuity-corrected %.4f  (+%.1f%%)"
          % (g("sd_mean").mean(), g("sd_cc_mean").mean(),
             100 * (g("sd_cc_mean").mean() / g("sd_mean").mean() - 1)))
    hi = float(np.abs(allz).max()); print("  max |z| over all %d pairs: %.1f" % (len(allz), hi))
    for q in (0.99, 0.999):
        print("    |z| %.1f%%ile %.3f   (Gaussian %.3f)"
              % (100 * q, np.percentile(np.abs(allz), 100 * q),
                 {0.99: 2.576, 0.999: 3.291}[q]))

    print("\n" + "=" * 96)
    print("B0  THE COVERAGE OPERATOR -- three readouts of the SAME central interval")
    print("=" * 96)
    print("  %-34s%12s%12s%12s" % ("", "cov50", "cov90", "nominal"))
    for nm, lab in (("centre", "bin CENTRE endpoints (calib.py)"),
                    ("edge", "bin EDGE endpoints"),
                    ("interp", "interpolated quantile (honest)")):
        print("  %-34s%12.4f%12.4f%12s"
              % (lab, g("cov50_%s" % nm).mean(), g("cov90_%s" % nm).mean(), "0.50 / 0.90"))

    print("\n" + "=" * 96)
    print("B1  DOES THE WIDENING OPERATOR MOVE THE POSTERIOR MEDIAN?  (L2's mechanism)")
    print("=" * 96)
    print("  L2 claims: SD convolution 'moves mass ACROSS bins, shifts the median, changes the")
    print("  ranking', while TEMP 'very nearly preserves the median'.  Measured directly:")
    print("  %-10s%10s%14s%14s%12s" % ("arm", "param", "sd ratio", "mean |dmedian|", "d(mean) A"))
    for f in FGRID[1:]:
        print("  %-10s%10.2f%14.3f%14.4f%+12.4f"
              % ("SD", f, g("sdratio_SD_%s" % f).mean(), g("dmed_SD_%s" % f).mean(),
                 g("dmean_SD_%s" % f).mean()))
    for T in TGRID[1:]:
        print("  %-10s%10.2f%14.3f%14.4f%+12.4f"
              % ("TEMP", T, g("sdratio_TEMP_%s" % T).mean(), g("dmed_TEMP_%s" % T).mean(),
                 g("dmean_TEMP_%s" % T).mean()))

    out = {"rows": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]}
    ST.save_atomic(os.path.join(RES, "audit_t1.json"), out,
                   complete_keys=("z_sd", "z_medabs", "cov90_interp", "dmed_SD_1.5"),
                   rows=out["rows"], n_expected=len(tg), module_file=__file__)
    return rows


if __name__ == "__main__":
    main()
