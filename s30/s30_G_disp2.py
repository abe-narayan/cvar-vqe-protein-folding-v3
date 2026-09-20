#!/usr/bin/env python
"""s30/s30_G_disp2.py -- LANE G, Q1 ADDENDUM: the adversarial checks on my own positive.

F-G1 fired.  Four things must be true before that is reportable, and each is a way the result
could be an artefact rather than the literature's claim:

  (1) DROP-TOP.  10 targets carry 67% of the whole-sample gain.  If the dispersion split is just
      those ten relabelled, it collapses when they are removed.  (Rank statistics already resist
      this -- rho is computed on ranks -- but the SPLIT is a mean and does not.)
  (2) CHAIN LENGTH.  rho(n, DISP) = +0.27.  Partial it out.
  (3) IS DIVERGENCE THE TAIL?  The coordinator's inference "concentrates on divergent pools =>
      tail intervention" needs divergence => tail.  Measure it directly.
  (4) MEDIATION.  DISP -> contraction -> d, with the artifact measured, not assumed.

Reads only `s30/results/s30_G_disp.json`.  Nothing is recomputed.
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
SRC = os.path.join(HERE, "results", "s30_G_disp.json")
OUT = os.path.join(HERE, "results", "s30_G_disp2.json")
SEED = 30_09_20


def rank(x):
    from scipy.stats import rankdata
    return rankdata(np.asarray(x, float))


def spear(x, y):
    a, b = rank(x) - rank(x).mean(), rank(y) - rank(y).mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else float("nan")


def partial_spear(x, y, z):
    """Spearman(x, y) with rank(z) partialled out of both."""
    rx, ry, rz = rank(x), rank(y), rank(z)
    rz = np.column_stack([np.ones_like(rz), rz])
    ex = rx - rz @ np.linalg.lstsq(rz, rx, rcond=None)[0]
    ey = ry - rz @ np.linalg.lstsq(rz, ry, rcond=None)[0]
    d = np.sqrt((ex * ex).sum() * (ey * ey).sum())
    return float((ex * ey).sum() / d) if d > 0 else float("nan")


def perm_null(fn, d, disp, folds, n=10000, seed=SEED):
    """Fold-preserving permutation null for any statistic fn(d, disp_permuted)."""
    rng = np.random.default_rng(seed)
    folds = np.asarray(folds)
    obs = fn(d, disp)
    idx = {f: np.where(folds == f)[0] for f in sorted(set(folds.tolist()))}
    draws = np.empty(n)
    for b in range(n):
        p = np.asarray(disp, float).copy()
        for f, ix in idx.items():
            p[ix] = np.asarray(disp, float)[rng.permutation(ix)]
        draws[b] = fn(d, p)
    sd = float(draws.std(ddof=1))
    return dict(observed=float(obs), null_mean=float(draws.mean()), null_sd=sd,
                mde=float(2.8016 * sd), ratio=float(obs / (2.8016 * sd)) if sd else float("nan"),
                p_two_sided=float((np.abs(draws - draws.mean()) >= abs(obs - draws.mean())).mean()))


def fold_boot_gap(d, folds, mask, n=4000, seed=SEED):
    rng = np.random.default_rng(seed)
    d = np.asarray(d, float); folds = np.asarray(folds); mask = np.asarray(mask, bool)
    eff = float(d[mask].mean() - d[~mask].mean())
    uf = sorted(set(folds.tolist()))
    draws = []
    for _ in range(n):
        pick = rng.choice(uf, size=len(uf), replace=True)
        a = np.concatenate([d[mask & (folds == f)] for f in pick])
        b = np.concatenate([d[~mask & (folds == f)] for f in pick])
        if a.size >= 2 and b.size >= 2:
            draws.append(a.mean() - b.mean())
    draws = np.asarray(draws)
    se = float(draws.std(ddof=1))
    per = []
    for f in uf:
        a, b = d[mask & (folds == f)], d[~mask & (folds == f)]
        per.append(float(a.mean() - b.mean()) if a.size and b.size else float("nan"))
    ok = [p for p in per if np.isfinite(p)]
    return dict(effect=eff, n_hi=int(mask.sum()), n_lo=int((~mask).sum()),
                mean_hi=float(d[mask].mean()), mean_lo=float(d[~mask].mean()),
                se_fold=se, mde=float(2.8016 * se),
                ratio=float(eff / (2.8016 * se)) if se else float("nan"),
                fold_ci=[float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
                excludes_zero=bool(np.percentile(draws, 2.5) > 0 or np.percentile(draws, 97.5) < 0),
                folds_same_sign=f"{sum(1 for p in ok if np.sign(p) == np.sign(eff))}/{len(ok)}",
                per_fold=[None if not np.isfinite(p) else round(float(p), 4) for p in per])


def main():
    S = json.load(open(SRC))
    rows = S["rows"]
    d = np.array([r["d"] for r in rows])
    folds = np.array([r["fold"] for r in rows], int)
    n = np.array([r["n"] for r in rows], float)
    f18 = np.array([r["fail18"] for r in rows], bool)
    disp = {t: np.array([r[t] for r in rows], float) for t in ("DISP_rmsd", "DISP_S", "DISP_rg")}
    contr = np.array([r["contraction"] for r in rows], float)
    rama = np.array([r["raw_rama_ok"] for r in rows], float)
    prod = np.array([r["proj"] for r in rows], float)
    poolmean = np.array([r["pool_mean_oracle"] for r in rows], float)

    out = dict(prereg="s30/PREREG_S30_G.md @ 78b65521", src=os.path.basename(SRC), n=len(rows))

    # --- (1) DROP-TOP: is the split just the ten winners relabelled?
    order = np.argsort(d)                       # most negative (biggest winners) first
    out["drop_top"] = {}
    for kdrop in (0, 5, 10, 20):
        keep = np.ones(len(d), bool); keep[order[:kdrop]] = False
        sub = {}
        for t, v in disp.items():
            m = v[keep] > np.median(v[keep])
            g = fold_boot_gap(d[keep], folds[keep], m)
            g["rho"] = spear(d[keep], v[keep])
            sub[t] = g
        sub["whole_sample_mean"] = float(d[keep].mean())
        out["drop_top"][f"drop_top{kdrop}"] = sub

    # --- winners' dispersion percentile: where do the top 10 sit?
    out["top10_winners"] = dict(
        pdbs=[rows[i]["pdb"] for i in order[:10]],
        d=[round(float(d[i]), 4) for i in order[:10]],
        disp_rmsd_pctile=[round(float((disp["DISP_rmsd"] < disp["DISP_rmsd"][i]).mean()), 3)
                          for i in order[:10]],
        median_pctile=float(np.median([(disp["DISP_rmsd"] < disp["DISP_rmsd"][i]).mean()
                                       for i in order[:10]])))

    # --- (2) CHAIN LENGTH partialled out
    out["partial_n"] = {}
    for t, v in disp.items():
        out["partial_n"][t] = dict(
            raw=spear(d, v), partial_n=partial_spear(d, v, n),
            null=perm_null(lambda a, b: partial_spear(a, b, n), d, v, folds, n=4000))
    out["partial_n"]["rho(n, d)"] = spear(n, d)
    out["partial_n"]["rho(n, DISP_rmsd)"] = spear(n, disp["DISP_rmsd"])
    # and the length-matched split: high/low dispersion WITHIN length tertiles
    ter = np.digitize(n, np.quantile(n, [1 / 3, 2 / 3]))
    m = np.zeros(len(d), bool)
    for q in (0, 1, 2):
        s = ter == q
        if s.sum() >= 4:
            m[s] = disp["DISP_rmsd"][s] > np.median(disp["DISP_rmsd"][s])
    out["length_matched_split"] = fold_boot_gap(d, folds, m)
    out["length_matched_split"]["note"] = "high vs low dispersion WITHIN chain-length tertiles"

    # --- (3) IS DIVERGENCE THE TAIL?
    thr = np.sort(poolmean)[-18] if np.isfinite(poolmean).all() else np.nan
    out["is_divergence_the_tail"] = dict(
        rho_disp_vs_production_rmsd=spear(disp["DISP_rmsd"], prod),
        rho_disp_vs_pool_mean_ORACLE=spear(disp["DISP_rmsd"], poolmean),
        FAIL18_dispersion_contrast=fold_boot_gap(disp["DISP_rmsd"], folds, f18),
        ORACLE_worst18_poolmean_dispersion_contrast=(
            fold_boot_gap(disp["DISP_rmsd"], folds, poolmean >= thr) if np.isfinite(thr) else None),
        overlap_hi_disp_and_FAIL18=dict(
            n_fail18_in_hi_disp=int((f18 & (disp["DISP_rmsd"] > np.median(disp["DISP_rmsd"]))).sum()),
            n_fail18=int(f18.sum()), expected_if_independent=9.0),
        note="ORACLE rows labelled ORACLE. FAIL18 is defined by the FILTER's in-band recall "
             "(S30-L23), so it is a descriptive stratum here.")

    # --- (4) MEDIATION: DISP -> artifact -> d
    out["mediation"] = dict(
        rho_disp_contraction=spear(disp["DISP_rmsd"], contr),
        rho_disp_rama=spear(disp["DISP_rmsd"], rama),
        rho_contraction_d=spear(contr, d),
        rho_disp_d=spear(disp["DISP_rmsd"], d),
        rho_disp_d_given_contraction=partial_spear(d, disp["DISP_rmsd"], contr),
        rho_contraction_d_given_disp=partial_spear(d, contr, disp["DISP_rmsd"]),
        note="if DISP acts ONLY through the artifact, partialling contraction kills rho(disp,d) "
             "and partialling disp leaves rho(contraction,d) -- that IS the literature's mechanism, "
             "not a confound.")

    # --- endpoint arithmetic: what is this worth, and does it touch the opening counterfactual?
    hi = disp["DISP_rmsd"] > np.median(disp["DISP_rmsd"])
    out["endpoint_arithmetic"] = dict(
        whole_sample_gain=float(d.mean()),
        gain_if_applied_only_to_hi_half=float(d[hi].sum() / len(d)),
        share_of_gain_in_hi_half=float(d[hi].sum() / d.sum()),
        production_mean=3.2105, pct_of_baseline=float(abs(d.mean()) / 3.2105 * 100),
        FAIL18_mean_d=float(d[f18].mean()), other108_mean_d=float(d[~f18].mean()),
        note="the opening counterfactual asks for -0.30 A from capping the worst 10 at 3.00 A")

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    print(f"wrote {OUT}\n")
    print("=== (1) DROP-TOP (does the split survive removing the biggest winners?) ===")
    for k, v in out["drop_top"].items():
        g = v["DISP_rmsd"]
        print(f"{k:10s} whole {v['whole_sample_mean']:+.4f} | hi-lo {g['effect']:+.4f} "
              f"[{g['fold_ci'][0]:+.4f},{g['fold_ci'][1]:+.4f}] {g['ratio']:+.2f}x MDE "
              f"folds {g['folds_same_sign']} | rho {g['rho']:+.3f}")
    w = out["top10_winners"]
    print(f"top-10 winners' dispersion percentile: median {w['median_pctile']:.3f}  {w['disp_rmsd_pctile']}")

    print("\n=== (2) CHAIN LENGTH ===")
    for t in ("DISP_rmsd", "DISP_S", "DISP_rg"):
        p = out["partial_n"][t]
        print(f"{t:10s} rho {p['raw']:+.3f} -> partial-n {p['partial_n']:+.3f} "
              f"({p['null']['ratio']:+.2f}x MDE, p2 {p['null']['p_two_sided']:.3f})")
    print(f"rho(n,d) {out['partial_n']['rho(n, d)']:+.3f}  rho(n,DISP) {out['partial_n']['rho(n, DISP_rmsd)']:+.3f}")
    lm = out["length_matched_split"]
    print(f"length-matched split: {lm['effect']:+.4f} [{lm['fold_ci'][0]:+.4f},{lm['fold_ci'][1]:+.4f}] "
          f"{lm['ratio']:+.2f}x MDE folds {lm['folds_same_sign']}")

    print("\n=== (3) IS DIVERGENCE THE TAIL? ===")
    t3 = out["is_divergence_the_tail"]
    print(f"rho(DISP, production RMSD) {t3['rho_disp_vs_production_rmsd']:+.3f}   "
          f"rho(DISP, pool mean ORACLE) {t3['rho_disp_vs_pool_mean_ORACLE']:+.3f}")
    f = t3["FAIL18_dispersion_contrast"]
    print(f"FAIL18 dispersion vs the 108: {f['effect']:+.4f} [{f['fold_ci'][0]:+.4f},{f['fold_ci'][1]:+.4f}] "
          f"{f['ratio']:+.2f}x MDE")
    o = t3["ORACLE_worst18_poolmean_dispersion_contrast"]
    if o:
        print(f"ORACLE worst18-by-pool-mean dispersion: {o['effect']:+.4f} "
              f"[{o['fold_ci'][0]:+.4f},{o['fold_ci'][1]:+.4f}] {o['ratio']:+.2f}x MDE")
    print(f"FAIL18 in the high-dispersion half: {t3['overlap_hi_disp_and_FAIL18']['n_fail18_in_hi_disp']}"
          f"/18 (9 expected if independent)")

    print("\n=== (4) MEDIATION ===")
    print(json.dumps({k: (round(v, 3) if isinstance(v, float) else v)
                      for k, v in out["mediation"].items() if k != "note"}))
    print("\n=== ENDPOINT ===")
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v)
                      for k, v in out["endpoint_arithmetic"].items() if k != "note"}))


if __name__ == "__main__":
    main()
