"""S25 / PHYSICS LANE -- STATISTICS FOR THE 7-CONFIGURATION SUITE. Pre-reg: `s25/PREREG_PHYS.md`.

    python s25/phys_analyse.py

Reads `s25/results/suite_cells/*.json` (written by `phys_suite.py`), asserts the two soundness
falsifiers, and emits every statistic the brief requires: per configuration and per pairwise
contrast, mean / median / SD / quartiles / best / worst target, paired delta, SE,
MDE = 2.8016 x SE computed PER COMPARISON, effect/MDE, W/L, iid CI beside fold-clustered CI,
worst-target degradation, and fold-wise plus length-stratified breakdowns.

Everything goes through `s24/stats_lib.py` verbatim. Nothing is hand-rolled -- the whole point
of that module is that two lanes' numbers can be laid side by side without an operator audit,
and its verdict rule has already been corrected once (it labelled 0.39x MDE "MEASURED").

THE TWO GATES THAT RUN FIRST, AND STOP EVERYTHING IF THEY FAIL
==============================================================
F1  Configuration 3 under the classical top-75 control must return the incumbent 3.0483 A.
    `zrank` is monotone, so this is an EXACTNESS check, not an approximation.
F2  Subset-hood (the s24 set-equality theorem) on every cell, and |rmsd_vqe - rmsd_topm| ~ 0.
    A non-zero eps is a bug, never a quantum effect.

BEST-OF-SEVEN IS AN ORDER STATISTIC
===================================
"the best configuration" chosen by looking at seven dev means is a minimum over seven
correlated draws. It is priced against `ST.best_of_k_null` with `share_accounted` and an
effective k that accounts for the correlation between configurations, and it is never quoted
as though it were a measured selection rule.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402

CELLS = os.path.join(P.RESULTS, "suite_cells")
INCUMBENT = 3.0483380938795324             # s21/results/poolgap.json avg_75, the endpoint


def load() -> List[Dict]:
    rows = []
    for f in sorted(glob.glob(os.path.join(CELLS, "*.json"))):
        rows.extend(json.load(open(f))["rows"])
    return rows


def grid(rows, norm, arm):
    """(126,) x 7 matrix of the chosen arm, targets in pinned order, configs 1..7."""
    pdbs = sorted({r["pdb"] for r in rows})
    idx = {p: i for i, p in enumerate(pdbs)}
    M = np.full((len(pdbs), 7), np.nan)
    for r in rows:
        if r["norm"] == norm:
            M[idx[r["pdb"]], r["config"] - 1] = r[arm]
    return pdbs, M


def desc(v):
    v = np.asarray(v, float)
    q = np.percentile(v, [25, 50, 75])
    return dict(mean=float(v.mean()), median=float(q[1]), sd=float(v.std(ddof=1)),
                q1=float(q[0]), q3=float(q[2]), min=float(v.min()), max=float(v.max()))


def main() -> int:
    rows = load()
    pdbs = sorted({r["pdb"] for r in rows})
    n = len(pdbs)
    meta = {r["pdb"]: (r["fold"], r["n"]) for r in rows}
    folds = np.array([meta[p][0] for p in pdbs], int)
    lens = np.array([meta[p][1] for p in pdbs], int)
    print(f"S25 PHYSICS LANE -- 7-CONFIGURATION COMPARISON SUITE")
    print(f"n = {n} targets x {len(P.CONFIGS)} configurations x {len(('rank','moment'))} "
          f"normalisations = {len(rows)} cells")
    print(f"basis = point_cloud   metric = full-chain CA-RMSD   selector = CVaR-VQE "
          f"(alpha={P.ALPHA}, T={P.TEMP}, {P.LAYERS} layers, {P.ITERS} iters, seed {P.SEED})")
    assert len(rows) == n * 14, f"incomplete: {len(rows)} rows for {n} targets"

    out = dict(n=n, alpha=P.ALPHA, T=P.TEMP, layers=P.LAYERS, iters=P.ITERS, seed=P.SEED,
               m_prod=P.M_PROD, k=P.K, incumbent=INCUMBENT, basis="point_cloud")

    # ============================================================ F1, F2 -- the gates
    print("\n" + "=" * 92)
    print("GATES")
    print("=" * 92)
    _, T75 = grid(rows, "rank", "rmsd_top75")
    anchor = float(np.nanmean(T75[:, 2]))
    f1 = abs(anchor - INCUMBENT) < 5e-5
    print(f"  F1 ANCHOR    C3 (Distogram) classical top-75 = {anchor:.6f} A   "
          f"incumbent {INCUMBENT:.6f}   |diff| {abs(anchor-INCUMBENT):.2e}   "
          f"{'PASS' if f1 else 'FAIL'}")
    gp = all(r["gate_pass"] for r in rows)
    eps = np.array([r["eps_vqe_minus_topm"] for r in rows], float)
    ge = float(np.mean([r["gate_equality"] for r in rows]))
    f2 = gp and float(np.abs(eps).max()) < 1e-9
    print(f"  F2 THEOREM   subset-hood {sum(r['gate_pass'] for r in rows)}/{len(rows)} cells   "
          f"max |rmsd_vqe - rmsd_topm| = {np.abs(eps).max():.2e}   {'PASS' if f2 else 'FAIL'}")
    print(f"               equality (reported, not asserted): {ge:.4f} of cells; "
          f"mean holes {np.mean([r['gate_n_holes'] for r in rows]):.2f}")
    out["gates"] = dict(f1_anchor=anchor, f1_pass=bool(f1), f2_pass=bool(f2),
                        max_abs_eps=float(np.abs(eps).max()), equality_rate=ge)
    if not (f1 and f2):
        print("\n  GATE FAILURE -- nothing downstream is read.")
        return 1

    # ============================================================ the headline table
    for norm in ("rank", "moment"):
        _, V = grid(rows, norm, "rmsd_vqe")
        _, C75 = grid(rows, norm, "rmsd_top75")
        _, CM = grid(rows, norm, "rmsd_topm")
        _, PM = grid(rows, norm, "rmsd_perm75")
        _, MQ = grid(rows, norm, "q_m")
        _, EN = grid(rows, norm, "q_entropy_bits")
        _, ES = grid(rows, norm, "q_ess")
        rnd = np.array([[r["rand_mean"] for r in rows
                         if r["pdb"] == p and r["norm"] == norm][0] for p in pdbs])
        rndb = np.array([[r["rand_best"] for r in rows
                          if r["pdb"] == p and r["norm"] == norm][0] for p in pdbs])
        tag = "PRIMARY (rank standardisation)" if norm == "rank" else \
              "DECLARED SECONDARY (raw moment z)"
        print("\n" + "=" * 92)
        print(f"THE SEVEN CONFIGURATIONS -- {tag}")
        print("=" * 92)
        print(f"  {'#':>2} {'configuration':<24} {'VQE':>8} {'median':>8} {'SD':>7} "
              f"{'Q1':>7} {'Q3':>7} {'best':>7} {'worst':>7} {'m':>5} {'ent':>6}")
        cfg_stats = {}
        for cid, name, subset in P.CONFIGS:
            v = V[:, cid - 1]
            d = desc(v)
            b, w = pdbs[int(np.argmin(v))], pdbs[int(np.argmax(v))]
            print(f"  {cid:>2} {name:<24} {d['mean']:>8.4f} {d['median']:>8.4f} "
                  f"{d['sd']:>7.4f} {d['q1']:>7.4f} {d['q3']:>7.4f} "
                  f"{d['min']:>7.3f} {d['max']:>7.3f} {MQ[:,cid-1].mean():>5.1f} "
                  f"{EN[:,cid-1].mean():>6.2f}")
            cfg_stats[cid] = dict(name=name, subset="+".join(subset), **d,
                                  best_target=b, worst_target=w,
                                  mean_m=float(MQ[:, cid - 1].mean()),
                                  mean_entropy_bits=float(EN[:, cid - 1].mean()),
                                  mean_ess=float(ES[:, cid - 1].mean()))
        print(f"\n  ARMS, per configuration (all POINT-CLOUD, all the SAME readout)")
        print(f"  {'#':>2} {'configuration':<24} {'VQE':>8} {'top-75':>8} {'top-m':>8} "
              f"{'perm75':>8} {'m-ladder':>9}")
        for cid, name, _ in P.CONFIGS:
            pm = PM[:, cid - 1]
            pms = f"{np.nanmean(pm):>8.4f}" if np.isfinite(pm).any() else f"{'--':>8}"
            print(f"  {cid:>2} {name:<24} {V[:,cid-1].mean():>8.4f} "
                  f"{C75[:,cid-1].mean():>8.4f} {CM[:,cid-1].mean():>8.4f} {pms} "
                  f"{(CM[:,cid-1]-C75[:,cid-1]).mean():>+9.4f}")
            cfg_stats[cid].update(mean_top75=float(C75[:, cid - 1].mean()),
                                  mean_topm=float(CM[:, cid - 1].mean()),
                                  mean_perm75=float(np.nanmean(pm)) if np.isfinite(pm).any()
                                  else None,
                                  m_ladder=float((CM[:, cid - 1] - C75[:, cid - 1]).mean()))
        print(f"  {'--':>2} {'random 75-subset (null)':<24} {rnd.mean():>8.4f}"
              f"   [best-of-16 {rndb.mean():.4f}]")
        out[f"configs_{norm}"] = cfg_stats
        out[f"random_null_{norm}"] = dict(mean=float(rnd.mean()),
                                          best_of_16=float(rndb.mean()))
        if norm == "rank":
            Vr, C75r, CMr, PMr, MQr, rndr = V, C75, CM, PM, MQ, rnd

    # ============================================================ contrasts
    print("\n" + "=" * 92)
    print("EVERY PAIRWISE CONTRAST, PAIRED PER TARGET -- PRIMARY (rank), VQE ARM")
    print("  d = a - b; NEGATIVE means the first configuration is BETTER")
    print("=" * 92)
    pair = {}
    for i in range(7):
        for j in range(i + 1, 7):
            a, b = P.CONFIGS[i], P.CONFIGS[j]
            lab = f"C{a[0]} {a[1]}  vs  C{b[0]} {b[1]}"
            r = ST.compare(Vr[:, i], Vr[:, j], folds, names=pdbs, label=lab)
            pair[f"C{a[0]}_vs_C{b[0]}"] = r
            print(f"\n{ST.fmt(r)}")
    out["pairwise_rank_vqe"] = pair

    print("\n" + "=" * 92)
    print("EACH CONFIGURATION AGAINST THE INCUMBENT (C3 classical top-75, 3.0483 A)")
    print("=" * 92)
    base = T75[:, 2]
    vs_inc = {}
    for cid, name, _ in P.CONFIGS:
        r = ST.compare(Vr[:, cid - 1], base, folds, names=pdbs,
                       label=f"C{cid} {name} (VQE)  vs  INCUMBENT")
        vs_inc[cid] = r
        print(f"\n{ST.fmt(r)}")
    out["vs_incumbent"] = vs_inc

    # ============================================================ VQE vs its controls
    print("\n" + "=" * 92)
    print("CVaR-VQE AGAINST ITS MATCHED CONTROLS, PER CONFIGURATION (rank)")
    print("=" * 92)
    ctrl = {}
    for cid, name, _ in P.CONFIGS:
        i = cid - 1
        c = {}
        c["vs_top75"] = ST.compare(Vr[:, i], C75r[:, i], folds, names=pdbs,
                                   label=f"C{cid} VQE vs classical top-75")
        c["vs_topm"] = ST.compare(Vr[:, i], CMr[:, i], folds, names=pdbs,
                                  label=f"C{cid} VQE vs classical top-m (size-matched)")
        c["vs_random"] = ST.compare(Vr[:, i], rndr, folds, names=pdbs,
                                    label=f"C{cid} VQE vs random 75-subset")
        if np.isfinite(PMr[:, i]).all():
            c["vs_perm"] = ST.compare(C75r[:, i], PMr[:, i], folds, names=pdbs,
                                      label=f"C{cid} top-75 vs RANK-PERMUTED physics channel")
        ctrl[cid] = c
        print(f"\n  --- C{cid} {name}")
        for k2, v in c.items():
            print(f"    {k2:<12} {v['effect']:+.4f}  SE {v['se']:.4f}  "
                  f"{v['effect_over_mde']:+.2f}x MDE  fold CI "
                  f"[{v['ci95_fold'][0]:+.4f},{v['ci95_fold'][1]:+.4f}]  "
                  f"{v['n_better']}W/{v['n_worse']}L  {v['verdict']}")
    out["controls_rank"] = ctrl

    # ============================================================ normalisation arm
    print("\n" + "=" * 92)
    print("F4 -- THE NORMALISATION FORK: rank vs raw moment, PAIRED, VQE arm")
    print("=" * 92)
    _, Vm = grid(rows, "moment", "rmsd_vqe")
    nrm = {}
    for cid, name, _ in P.CONFIGS:
        r = ST.compare(Vm[:, cid - 1], Vr[:, cid - 1], folds, names=pdbs,
                       label=f"C{cid} {name}: MOMENT vs RANK")
        nrm[cid] = r
        print(f"  C{cid} {name:<24} moment {Vm[:,cid-1].mean():.4f}  rank "
              f"{Vr[:,cid-1].mean():.4f}  d {r['effect']:+.4f}  "
              f"{r['effect_over_mde']:+.2f}x MDE  {r['n_better']}W/{r['n_worse']}L  "
              f"{r['verdict']}")
    out["normalisation_fork"] = nrm

    # ============================================================ breakdowns
    print("\n" + "=" * 92)
    print("FOLD-WISE AND LENGTH-STRATIFIED BREAKDOWN (rank, VQE arm)")
    print("=" * 92)
    print(f"  {'#':>2} {'configuration':<24} " +
          " ".join(f"{'f%d' % f:>7}" for f in sorted(set(folds.tolist()))))
    fw = {}
    for cid, name, _ in P.CONFIGS:
        v = Vr[:, cid - 1]
        per = {int(f): float(v[folds == f].mean()) for f in sorted(set(folds.tolist()))}
        fw[cid] = per
        print(f"  {cid:>2} {name:<24} " + " ".join(f"{per[f]:>7.4f}" for f in sorted(per)))
    edges = [0, 11, 13, 15, 99]
    bands = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    print(f"\n  {'#':>2} {'configuration':<24} " +
          " ".join(f"{'n<=%d' % b[1] if b[1] < 99 else 'n>%d' % b[0]:>9}" for b in bands))
    lw = {}
    for cid, name, _ in P.CONFIGS:
        v = Vr[:, cid - 1]
        per = {f"{lo}-{hi}": float(v[(lens > lo) & (lens <= hi)].mean()) for lo, hi in bands}
        lw[cid] = per
        print(f"  {cid:>2} {name:<24} " +
              " ".join(f"{per[f'{lo}-{hi}']:>9.4f}" for lo, hi in bands))
    print("  band sizes: " + " ".join(f"{lo}-{hi}:{int(((lens>lo)&(lens<=hi)).sum())}"
                                      for lo, hi in bands))
    out["fold_breakdown"] = fw
    out["length_breakdown"] = lw

    # ============================================================ best-of-7 order statistic
    print("\n" + "=" * 92)
    print("'THE BEST OF THE SEVEN' IS AN ORDER STATISTIC, AND IS PRICED AS ONE")
    print("=" * 92)
    means = Vr.mean(0)
    obs_best = float(means.min())
    best_cid = int(np.argmin(means)) + 1
    # per-target min over the seven: a fully leaked ORACLE selection rule
    per_t_min = Vr.min(1)
    Cm = np.corrcoef(Vr.T)
    off = Cm[np.triu_indices(7, 1)]
    k_eff = float(7.0 / (1.0 + (7 - 1) * float(off.mean())))
    bok = ST.best_of_k_null(Vr.ravel(), 7, minimise=True)
    print(f"  best configuration by dev mean: C{best_cid} at {obs_best:.4f} A")
    print(f"  mean pairwise correlation between the seven arms: {off.mean():.4f}  "
          f"-> k_eff = {k_eff:.2f} independent configurations, not 7")
    print(f"  ORACLE per-target min over the seven: {per_t_min.mean():.4f} A "
          f"({per_t_min.mean()-obs_best:+.4f} vs the best single configuration)")
    nb = ST.best_of_k_accounted(per_t_min.mean(), obs_best, Vr.ravel(), 7)
    print(f"  best-of-7 null on the pooled configuration distribution: gain "
          f"{nb['null_gain']:+.4f} vs observed {nb['observed_gain']:+.4f}  "
          f"share_accounted {nb['share_accounted']:.3f}  residual {nb['residual_gain']:+.4f}")
    print(f"  -> the per-target min is ORACLE and is NOT a selection rule.")
    out["best_of_seven"] = dict(best_config=best_cid, best_mean=obs_best,
                                mean_pairwise_corr=float(off.mean()), k_eff=k_eff,
                                oracle_per_target_min=float(per_t_min.mean()),
                                accounted=nb)

    ST.save_atomic(os.path.join(P.RESULTS, "phys_suite.json"), out,
                   complete_keys=("pdb", "config", "norm", "rmsd_vqe", "rmsd_top75",
                                  "rmsd_topm", "gate_pass"),
                   rows=rows, n_expected=n * 14, module_file=__file__)
    print(f"\nwrote s25/results/phys_suite.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
