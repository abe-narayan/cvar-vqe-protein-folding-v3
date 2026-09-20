#!/usr/bin/env python
"""s30/s30_G_disp.py -- LANE G, Q1: where does the k=30 AMBER-relax gain live?

Pre-registered in `s30/PREREG_S30_G.md` (commit 78b65521).  NOTHING IS RECOMPUTED: the relax
arm is read from `s16/results/repair_A.json`, whose `k30_full` setting reproduces the project
memory's -0.022 [-0.036,-0.009] exactly (`s16/results/repair_report.json`).

The literature claim under test (lane L: asserted, never measured): "averaging artifacts become
more pronounced when members of the ensemble are more divergent."

    python s30/s30_G_disp.py            # writes s30/results/s30_G_disp.json
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
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPAIR = os.path.join(ROOT, "s16", "results", "repair_A.json")
REPORT = os.path.join(ROOT, "s16", "results", "repair_report.json")
XPOOL = os.path.join(RESULTS, "s30_X_typicalgood.json")
FSCORE = os.path.join(RESULTS, "s30_F_score.json")
OUT = os.path.join(RESULTS, "s30_G_disp.json")
SEED = 30_09_20


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return float("nan")
    from scipy.stats import rankdata
    a, b = rankdata(x[ok]), rankdata(y[ok])
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else float("nan")


def mean_ci(vals, folds, label, n_boot=4000):
    from s24 import stats_lib as ST
    v = np.asarray(vals, float)
    ok = np.isfinite(v)
    o = ST.compare(v[ok], np.zeros(int(ok.sum())), folds=list(np.asarray(folds)[ok]),
                   label=label, n_boot=n_boot)
    ci = o["ci95_fold"]
    return dict(mean=float(o["effect"]), median=float(o["median_effect"]), se=float(o["se"]),
                mde=float(o["mde"]),
                ratio=float(o["effect"] / o["mde"]) if o["mde"] else float("nan"),
                fold_ci=[float(ci[0]), float(ci[1])],
                folds_same_sign=o.get("folds_same_sign"), n=int(ok.sum()),
                excludes_zero=bool(ci[0] > 0 or ci[1] < 0))


def two_group(v, folds, mask, label):
    """Unpaired contrast mean(v[mask]) - mean(v[~mask]) with a fold-clustered CI, via the
    project's own comparator on the *stacked* samples (it does not require pairing)."""
    from s24 import stats_lib as ST
    v = np.asarray(v, float); folds = np.asarray(folds); mask = np.asarray(mask, bool)
    a, b = v[mask], v[~mask]
    fa, fb = folds[mask], folds[~mask]
    rng = np.random.default_rng(SEED)
    eff = float(a.mean() - b.mean())
    # fold-clustered bootstrap: resample FOLDS with replacement, recompute on the union
    uf = sorted(set(folds.tolist()))
    draws = []
    for _ in range(4000):
        pick = rng.choice(uf, size=len(uf), replace=True)
        va = np.concatenate([a[fa == f] for f in pick]) if any((fa == f).any() for f in pick) else np.array([])
        vb = np.concatenate([b[fb == f] for f in pick]) if any((fb == f).any() for f in pick) else np.array([])
        if va.size < 2 or vb.size < 2:
            continue
        draws.append(va.mean() - vb.mean())
    draws = np.asarray(draws, float)
    lo, hi = (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))) if draws.size > 50 else (float("nan"),) * 2
    se_iid = float(np.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size))
    se_fold = float(draws.std(ddof=1)) if draws.size > 50 else float("nan")
    se = max(se_iid, se_fold) if np.isfinite(se_fold) else se_iid
    mde = 2.8016 * se
    same = []
    for f in uf:
        aa, bb = a[fa == f], b[fb == f]
        same.append(float(aa.mean() - bb.mean()) if aa.size >= 1 and bb.size >= 1 else float("nan"))
    ss = [s for s in same if np.isfinite(s)]
    return dict(label=label, effect=eff, n_hi=int(a.size), n_lo=int(b.size),
                mean_hi=float(a.mean()), mean_lo=float(b.mean()),
                median_hi=float(np.median(a)), median_lo=float(np.median(b)),
                se_iid=se_iid, se_fold=se_fold, mde=float(mde),
                ratio=float(eff / mde) if mde else float("nan"),
                fold_ci=[lo, hi], excludes_zero=bool(np.isfinite(lo) and (lo > 0 or hi < 0)),
                per_fold=[None if not np.isfinite(s) else float(s) for s in same],
                folds_same_sign=f"{sum(1 for s in ss if np.sign(s) == np.sign(eff))}/{len(ss)}")


def perm_null_rho(d, disp, folds, n=10000):
    """Fold-PRESERVING permutation null for Spearman(d, disp): the dispersion labels are permuted
    WITHIN fold, so fold-level structure in either variable cannot manufacture a correlation."""
    rng = np.random.default_rng(SEED)
    folds = np.asarray(folds)
    obs = spearman(d, disp)
    draws = np.empty(n)
    idx_by_fold = {f: np.where(folds == f)[0] for f in sorted(set(folds.tolist()))}
    dd = np.asarray(disp, float)
    for b in range(n):
        p = dd.copy()
        for f, ix in idx_by_fold.items():
            p[ix] = dd[rng.permutation(ix)]
        draws[b] = spearman(d, p)
    sd = float(draws.std(ddof=1))
    return dict(rho=float(obs), null_mean=float(draws.mean()), null_sd=sd,
                mde_rho=float(2.8016 * sd), ratio=float(obs / (2.8016 * sd)) if sd else float("nan"),
                p_two_sided=float((np.abs(draws) >= abs(obs)).mean()),
                p_one_sided_neg=float((draws <= obs).mean()))


def uniform_effect_null(d, disp, folds, n=10000):
    """Contract rule 5.  Under a UNIFORM effect the per-target d_t are exchangeable with respect to
    dispersion, so permute the dispersion labels (within fold) and read where the observed
    high-minus-low gap sits in that null.  This is the valid test; a raw drop-top threshold is not."""
    rng = np.random.default_rng(SEED + 1)
    folds = np.asarray(folds); d = np.asarray(d, float); disp = np.asarray(disp, float)
    med = np.median(disp)
    hi = disp > med
    obs = float(d[hi].mean() - d[~hi].mean())
    idx_by_fold = {f: np.where(folds == f)[0] for f in sorted(set(folds.tolist()))}
    draws = np.empty(n)
    for b in range(n):
        p = disp.copy()
        for f, ix in idx_by_fold.items():
            p[ix] = disp[rng.permutation(ix)]
        h = p > np.median(p)
        draws[b] = d[h].mean() - d[~h].mean()
    return dict(observed_gap=obs, null_mean=float(draws.mean()), null_sd=float(draws.std(ddof=1)),
                pctile=float((draws <= obs).mean()),
                null_p2_5=float(np.percentile(draws, 2.5)),
                null_p97_5=float(np.percentile(draws, 97.5)))


def main():
    from s12 import instrument as I

    assert os.path.exists(REPAIR) and os.path.exists(REPORT), "artefact missing"
    A = json.load(open(REPAIR))
    R = json.load(open(REPORT))
    arm = R["settings"]["k30_full"]["vs_proj_ungated"]
    print(f"arm k30_full vs proj (ungated): mean {arm['mean_diff']:+.4f}  se {arm['se']:.4f}  "
          f"n {arm['n']}  W/L {arm['n_better']}/{arm['n_worse']}")

    xs = {r["pdb"]: r for r in json.load(open(XPOOL))["pool"]} if os.path.exists(XPOOL) else {}
    fs = {r["pdb"]: r for r in json.load(open(FSCORE))["per_target"]} if os.path.exists(FSCORE) else {}

    rows = []
    for t in A["per_target"]:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        pool = I.pool_idx(u)
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        W75 = np.asarray(u["W"], float)[pool[sub]]
        P = I.pairwise_rmsd(W75)
        iu = np.triu_indices(len(W75), 1)
        disp_rmsd = float(P[iu].mean())
        C, b = I.coordinate_average(W75, P)
        Wc = I.superpose_batch(W75, W75[b])
        S_spread = float(np.sqrt((((Wc - C[None]) ** 2).sum(-1)).mean()))
        rg = np.sqrt((((W75 - W75.mean(1, keepdims=True)) ** 2).sum(-1)).mean(1))
        g = t["raw_avg_geom"]
        rows.append(dict(
            pdb=pdb, fold=int(t["fold"]), n=int(t["n"]),
            fail18=pdb in I.FAIL18,
            proj=float(t["proj"]["rmsd"]),
            k30=float(t["f0_k30_full"]["rmsd"]),
            d=float(t["f0_k30_full"]["rmsd"] - t["proj"]["rmsd"]),
            converged=bool(t["f0_k30_full"]["converged"]),
            raw_avg=float(t["raw_avg_rmsd"]),
            raw_rama_ok=float(t["raw_avg_rama_ok"]),
            contraction=float(1.0 - g["ca_ca"] / 3.80),          # CA-CA ideal 3.80 A
            moved_ca=float(t["f0_k30_full"]["moved_ca"]),
            restraint_rmsd=float(t["f0_k30_full"]["restraint_rmsd"]),
            DISP_rmsd=disp_rmsd, DISP_S=S_spread, DISP_rg=float(rg.std(ddof=1)),
            X_S_spread=float(xs[pdb]["S_spread"]) if pdb in xs else float("nan"),
            F_top75_rg_sd=float(fs[pdb]["F4_top75_rg_sd"]) if pdb in fs else float("nan"),
            pool_mean_oracle=float(fs[pdb]["pool_mean"]) if pdb in fs else float("nan"),
            pool_best_oracle=float(fs[pdb]["pool_best"]) if pdb in fs else float("nan"),
        ))
        if len(rows) % 25 == 0:
            print(f"  {len(rows)}/126")

    d = np.array([r["d"] for r in rows])
    folds = np.array([r["fold"] for r in rows], int)
    out = dict(prereg="s30/PREREG_S30_G.md @ 78b65521",
               arm="s16/results/repair_A.json :: f0_k30_full vs proj (k=30, steps=0)",
               arm_reproduced=dict(mean=float(d.mean()), se_iid=float(d.std(ddof=1) / np.sqrt(len(d))),
                                   report_mean=float(arm["mean_diff"]), report_se=float(arm["se"])),
               n=len(rows), rows=rows)

    out["whole_sample"] = mean_ci(d, folds, "G.k30_vs_proj")
    out["POWER_registered"] = dict(
        sd_d=float(d.std(ddof=1)),
        mde_half_split=float(2.8016 * d.std(ddof=1) * np.sqrt(2 / 63)),
        mde_fail18=float(2.8016 * d.std(ddof=1) * np.sqrt(1 / 18 + 1 / 108)),
        note="registered BEFORE the split: a half-split cannot reach its own MDE unless "
             "MORE than 100% of a 0.022 A effect sits in one half")

    out["cross_check_dispersion"] = {
        "DISP_S_vs_laneX_S_spread": spearman([r["DISP_S"] for r in rows], [r["X_S_spread"] for r in rows]),
        "DISP_rg_vs_laneF_top75_rg_sd": spearman([r["DISP_rg"] for r in rows], [r["F_top75_rg_sd"] for r in rows]),
        "DISP_rmsd_vs_DISP_S": spearman([r["DISP_rmsd"] for r in rows], [r["DISP_S"] for r in rows]),
    }

    out["F_G1_continuous"] = {}
    out["F_G1b_split"] = {}
    out["uniform_effect_null"] = {}
    for tag in ("DISP_rmsd", "DISP_S", "DISP_rg"):
        v = np.array([r[tag] for r in rows], float)
        out["F_G1_continuous"][tag] = perm_null_rho(d, v, folds)
        out["F_G1b_split"][tag] = two_group(d, folds, v > np.median(v), f"G.split.{tag}")
        out["uniform_effect_null"][tag] = uniform_effect_null(d, v, folds)

    # --- strata
    f18 = np.array([r["fail18"] for r in rows], bool)
    out["F_G1c_strata"] = dict(FAIL18_minus_108=two_group(d, folds, f18, "G.FAIL18"))
    pm = np.array([r["pool_mean_oracle"] for r in rows], float)
    if np.isfinite(pm).all():
        thr = np.sort(pm)[-18]
        out["F_G1c_strata"]["ORACLE_worst18_by_pool_mean_minus_rest"] = two_group(
            d, folds, pm >= thr, "G.ORACLE.worst18_poolmean")
    out["F_G1c_strata"]["note"] = ("FAIL18 is defined by the FILTER's in-band recall (S30-L23), not "
                                   "by d_t; reported as a descriptive split. The worst18-by-pool-mean "
                                   "row is ORACLE and labelled ORACLE. fold 0 has no FAIL18 target.")

    # --- concentration verdict as ONE object (mean / median / win-rate + the uniform null)
    out["concentration"] = dict(
        mean=float(d.mean()), median=float(np.median(d)),
        sd=float(d.std(ddof=1)), win_rate=float((d < 0).mean()),
        n_better=int((d < 0).sum()), n_worse=int((d > 0).sum()),
        drop_top10=float(np.sort(d)[10:].mean()), drop_top20=float(np.sort(d)[20:].mean()),
        top10_share=float(np.sort(d)[:10].sum() / d.sum()) if d.sum() else float("nan"),
        note="d_t is LINEAR in RMSD so the mean is the primary; median and W/L travel beside it. "
             "The drop-top numbers are DESCRIPTIVE -- the valid concentration test is the "
             "uniform-effect null above (contract rule 5).")

    # --- mechanism sub-check: does divergence predict the ARTIFACT at all?
    out["mechanism"] = {
        "rho(DISP_rmsd, contraction)": spearman([r["DISP_rmsd"] for r in rows], [r["contraction"] for r in rows]),
        "rho(DISP_rmsd, raw_avg_rama_ok)": spearman([r["DISP_rmsd"] for r in rows], [r["raw_rama_ok"] for r in rows]),
        "rho(DISP_rmsd, moved_ca)": spearman([r["DISP_rmsd"] for r in rows], [r["moved_ca"] for r in rows]),
        "rho(DISP_rmsd, raw_avg_rmsd)": spearman([r["DISP_rmsd"] for r in rows], [r["raw_avg"] for r in rows]),
        "rho(contraction, d)": spearman([r["contraction"] for r in rows], list(d)),
        "rho(moved_ca, d)": spearman([r["moved_ca"] for r in rows], list(d)),
        "rho(n, DISP_rmsd)": spearman([r["n"] for r in rows], [r["DISP_rmsd"] for r in rows]),
    }

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {OUT}")

    w = out["whole_sample"]
    print(f"\nWHOLE SAMPLE  d = {w['mean']:+.4f}  fold CI [{w['fold_ci'][0]:+.4f},{w['fold_ci'][1]:+.4f}]  "
          f"{w['ratio']:.2f}x MDE  folds {w['folds_same_sign']}")
    p = out["POWER_registered"]
    print(f"POWER         sd {p['sd_d']:.4f}   MDE half-split {p['mde_half_split']:.4f}   "
          f"MDE FAIL18 {p['mde_fail18']:.4f}")
    for tag in ("DISP_rmsd", "DISP_S", "DISP_rg"):
        c = out["F_G1_continuous"][tag]; s = out["F_G1b_split"][tag]; u = out["uniform_effect_null"][tag]
        print(f"\n{tag:10s} rho(d,disp) {c['rho']:+.3f}  null sd {c['null_sd']:.3f}  "
              f"{c['ratio']:+.2f}x MDE  p2 {c['p_two_sided']:.3f}")
        print(f"{'':10s} split hi-lo {s['effect']:+.4f}  [{s['fold_ci'][0]:+.4f},{s['fold_ci'][1]:+.4f}]  "
              f"{s['ratio']:+.2f}x MDE  folds {s['folds_same_sign']}  (hi {s['mean_hi']:+.4f} lo {s['mean_lo']:+.4f})")
        print(f"{'':10s} uniform-effect null: gap at pctile {u['pctile']:.3f} of "
              f"[{u['null_p2_5']:+.4f},{u['null_p97_5']:+.4f}]")
    for k, v in out["F_G1c_strata"].items():
        if isinstance(v, dict):
            print(f"\n{k}: {v['effect']:+.4f} [{v['fold_ci'][0]:+.4f},{v['fold_ci'][1]:+.4f}] "
                  f"{v['ratio']:+.2f}x MDE  n {v['n_hi']}/{v['n_lo']}  folds {v['folds_same_sign']}")
    print("\nmechanism:", json.dumps({k: round(v, 3) for k, v in out["mechanism"].items()}))
    cc = out["concentration"]
    print(f"concentration: mean {cc['mean']:+.4f} median {cc['median']:+.4f} win {cc['win_rate']:.3f} "
          f"drop10 {cc['drop_top10']:+.4f} top10_share {cc['top10_share']:.3f}")


if __name__ == "__main__":
    main()
