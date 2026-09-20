#!/usr/bin/env python
"""s30/s30_R_agg.py -- LANE R: the verdict on F-R1, from the rows.

Pre-registered in `s30/PREREG_S30_R.md` (commit 7eabffee).  Nothing here selects a channel: the
whole table is reported and the max-over-channels is priced against a per-target sign-flip null.

    python s30/s30_R_agg.py [--rows GLOB] [--out PATH]
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
DELTA_LABELS = ("0-0.25", "0.25-0.5", "0.5-1", "1-2", "2-4", ">4")

# --- F-R1's registered bars (PREREG_S30_R.md section 4).  Not editable after the fact.
BAR_RHO, BAR_RHO_MARGIN = 0.25, 0.10
BAR_PREF, BAR_PREF_MARGIN = 0.65, 0.10


def load_rows(pattern=None):
    pats = [pattern] if pattern else [os.path.join(RESULTS, "s30_R_rows.s*.jsonl")]
    rows, seen = [], set()
    for p in pats:
        for f in sorted(glob.glob(p)):
            with open(f, encoding="utf-8") as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    r = json.loads(line)
                    if r.get("error") or r.get("skipped") or r["pdb"] in seen:
                        continue
                    seen.add(r["pdb"]); rows.append(r)
    return sorted(rows, key=lambda r: r["pdb"])


def mean_ci(vals, folds, label="", n_boot=4000):
    """Mean with a FOLD-CLUSTERED CI and the per-comparison MDE, via the project's own stats."""
    from s24 import stats_lib as ST
    v = np.asarray(vals, float)
    ok = np.isfinite(v)
    if ok.sum() < 5:
        return dict(mean=float("nan"), n=int(ok.sum()), note="NOT MEASURED (n < 5)")
    o = ST.compare(v[ok], np.zeros(int(ok.sum())), folds=list(np.asarray(folds)[ok]), label=label,
                   n_boot=n_boot)
    ci = o["ci95_fold"]
    return dict(mean=float(o["effect"]), median=float(o["median_effect"]), se=float(o["se"]),
                mde=float(o["mde"]), ratio=float(o["effect"] / o["mde"]) if o["mde"] else float("nan"),
                fold_ci=[float(ci[0]), float(ci[1])], power=float(o.get("power", float("nan"))),
                folds_same_sign=o.get("folds_same_sign"), n=int(ok.sum()),
                excludes_zero=bool(ci[0] > 0 or ci[1] < 0))


def paired_ci(a, b, folds, label=""):
    from s24 import stats_lib as ST
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5:
        return dict(effect=float("nan"), n=int(ok.sum()), note="NOT MEASURED (n < 5)")
    o = ST.compare(a[ok], b[ok], folds=list(np.asarray(folds)[ok]), label=label)
    ci = o["ci95_fold"]
    return dict(effect=float(o["effect"]), se=float(o["se"]), mde=float(o["mde"]),
                ratio=float(o["effect"] / o["mde"]) if o["mde"] else float("nan"),
                fold_ci=[float(ci[0]), float(ci[1])], power=float(o.get("power", float("nan"))),
                folds_same_sign=o.get("folds_same_sign"), n=int(ok.sum()),
                excludes_zero=bool(ci[0] > 0 or ci[1] < 0))


def signflip_pmax(M, n_boot=500, seed=20300930):
    """Max-over-channels per-target sign-flip null (S28-L34(e)).  M is (n_targets, n_channels).
    Returns the null's mean/p95 of the max |mean| and the p of the observed max."""
    rng = np.random.default_rng(seed)
    M = np.asarray(M, float)
    keep = np.isfinite(M).sum(0) >= 5          # drop channels that are absent almost everywhere
    M = M[:, keep]
    if M.shape[1] == 0:
        return dict(note="NOT MEASURED (no channel with >= 5 targets)")
    obs = np.nanmax(np.abs(np.nanmean(M, axis=0)))
    draws = np.empty(n_boot)
    for b in range(n_boot):
        s = rng.choice([-1.0, 1.0], size=(M.shape[0], 1))
        draws[b] = np.nanmax(np.abs(np.nanmean(M * s, axis=0)))
    return dict(observed_max=float(obs), null_mean=float(draws.mean()),
                null_p95=float(np.percentile(draws, 95)),
                p_max=float((draws >= obs).mean()))


def lfo_combination(rows, names, key="d_near", ctrl="d_pool"):
    """Prereg item 6.  A ridge fit on four folds, evaluated on the fifth: does a COMBINATION of
    channels prefer the near-native rung to production where no single channel does?  Priced
    against the pool-member control (the same weights applied to a real pool member) and a
    label-shuffle null.  Lower score = preferred, so a correct vote is d < 0."""
    rng = np.random.default_rng(3030)
    folds = np.array([r["fold"] for r in rows], int)
    X = np.array([[r["ch"].get(nm, {}).get(key) if r["ch"].get(nm, {}).get(key) is not None
                   else np.nan for nm in names] for r in rows], float)
    C = np.array([[r["ch"].get(nm, {}).get(ctrl) if r["ch"].get(nm, {}).get(ctrl) is not None
                   else np.nan for nm in names] for r in rows], float)
    keep = np.isfinite(X).all(1) & np.isfinite(C).all(1)
    X, C, folds = X[keep], C[keep], folds[keep]
    if X.shape[0] < 20:
        return dict(note="NOT MEASURED (too few complete rows)", n=int(X.shape[0]))
    out = {}
    for tag, shuffle in (("real", False), ("shuffled", True)):
        pn, pc = np.zeros(len(X)), np.zeros(len(X))
        for f in sorted(set(folds)):
            tr, te = folds != f, folds == f
            A = X[tr].copy()
            if shuffle:                      # destroy the near/production correspondence
                A = A * rng.choice([-1.0, 1.0], size=(A.shape[0], 1))
            # weights that make w.d_near as NEGATIVE as possible, unit-norm, ridge-stabilised
            G = A.T @ A + 1.0 * np.eye(A.shape[1])
            w = -np.linalg.solve(G, A.sum(0))
            nw = np.linalg.norm(w)
            w = w / nw if nw > 1e-12 else w
            pn[te] = X[te] @ w
            pc[te] = C[te] @ w
        out[tag] = dict(
            pref_near=float(np.mean((pn < 0) + 0.5 * (pn == 0))),
            pref_pool=float(np.mean((pc < 0) + 0.5 * (pc == 0))),
            margin=float(np.mean((pn < 0) + 0.5 * (pn == 0)) - np.mean((pc < 0) + 0.5 * (pc == 0))),
            folds=list(map(int, folds)), _pn=pn.tolist(), _pc=pc.tolist())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", default=None)
    ap.add_argument("--out", default=os.path.join(RESULTS, "s30_R_verdict.json"))
    a = ap.parse_args(argv)

    rows = load_rows(a.rows)
    if not rows:
        print("no rows"); return
    folds = [r["fold"] for r in rows]
    names = sorted({k for r in rows for k in r["ch"]})
    print(f"n = {len(rows)} targets, {len(names)} channels, folds "
          f"{ {f: folds.count(f) for f in sorted(set(folds))} }")

    o = dict(n=len(rows), channels=names, prereg="s30/PREREG_S30_R.md @ 7eabffee")

    # ---------------- instrument scope, printed BEFORE any verdict
    o["instrument"] = dict(
        rebuild_floor=mean_ci([r["rebuild_floor"] for r in rows], folds, "floor"),
        contamination=float(np.mean([r["contamination"] for r in rows])),
        pool_ref_channels=rows[0].get("pool_ref_channels"),
        n_near=float(np.mean([r["n_near"] for r in rows])),
        near_rmsd=float(np.mean([r["near_rmsd_mean"] for r in rows])),
        prod_rmsd=float(np.mean([r["rmsd_fixed"]["PROD"] for r in rows])),
        circbest_rmsd=float(np.mean([r["rmsd_fixed"]["circ_best"] for r in rows])),
        pool_best=float(np.mean([r["pool_best"] for r in rows])),
        anchorB_rmsd=float(np.mean([r["anchorB_rmsd"] for r in rows])),
        n_fallback=int(sum(1 for r in rows if r.get("near_fallback"))))

    # --- A2 realism flatness (can only weaken my own positives)
    o["A2_realism_flatness"] = {k: mean_ci([r["realism_flatness"].get(k) for r in rows], folds, k)
                                for k in ("RAMA", "EXVOL", "RG_DEV")}
    # --- SI audit: the size-matched twin must make the pure Rg functions constant
    def _si(k, new_key, old_key):
        v = [r["si_audit"][k].get(new_key, r["si_audit"][k].get(old_key))
             for r in rows if k in r.get("si_audit", {})]
        v = [x for x in v if x is not None]
        return float(np.mean(v)) if v else float("nan")
    o["SI_audit"] = {k: dict(rel_sd_SI=_si(k, "sd_SI_over_raw_scale", "rel_sd"),
                             rel_sd_raw=_si(k, "sd_raw_over_raw_scale", "raw_rel_sd"))
                     for k in ("RG_LAW", "RG_UNIV")}
    # --- D3 anchor confound, measured BEFORE the verdict is read
    ac = np.array([r["anchor_confound"] for r in rows], float)
    o["D3_anchor_confound"] = dict(mean=float(np.nanmean(ac)), median=float(np.nanmedian(ac)),
                                   p10=float(np.nanpercentile(ac, 10)),
                                   p90=float(np.nanpercentile(ac, 90)),
                                   share_above_0p5=float(np.nanmean(ac > 0.5)),
                                   caveat_triggers=bool(np.nanmedian(ac) > 0.5))
    # --- D1 locality decomposition
    o["D1_locality"] = {k: mean_ci([r["locality"][k] for r in rows], folds, k)
                        for k in ("r2_m_only", "r2_local", "r2_global", "d_local", "d_global")}

    # ---------------- the per-channel table
    tab = {}
    for nm in names:
        g = lambda k: [r["ch"].get(nm, {}).get(k, float("nan")) for r in rows]
        cell = dict(n=int(sum(1 for r in rows if nm in r["ch"])))
        for k in ("rho_A", "rho_A_part", "rho_B", "rho_B_part", "rho_ANCHOR", "rho_ANCHOR_part",
                  "rho_B_part_anchor", "rho_rg_A", "rho_rg_pool", "pref_near", "pref_pool",
                  "pref_oracle_circbest", "pref_native_rebuilt", "pctile_nat_in_A"):
            cell[k] = mean_ci(g(k), folds, f"{nm}.{k}")
        # F-R1 clause (i)'s registered contrast: ladder A's nativeness ordering vs ladder B's
        # distance-to-anchor ordering.
        cell["contrast_A_minus_ANCHOR"] = paired_ci(g("rho_A_part"), g("rho_ANCHOR_part"), folds,
                                                    f"{nm}.A-ANCHOR")
        # Reported BESIDE it, never substituted for it: the SAME ladder-B structures with the
        # SAME channel values, only the LABEL changed (distance to the native vs distance to the
        # anchor).  This is the sharpest form the control can take -- nothing varies but the
        # question being asked.
        cell["contrast_B_minus_ANCHOR"] = paired_ci(g("rho_B_part"), g("rho_ANCHOR_part"), folds,
                                                    f"{nm}.B-ANCHOR")
        cell["contrast_pref"] = paired_ci(g("pref_near"), g("pref_pool"), folds, f"{nm}.pref")
        # D2 resolution: concordance per |delta RMSD| bin, count-weighted per target
        cell["conc"] = {}
        for q, lab in enumerate(DELTA_LABELS):
            v = [r["ch"].get(nm, {}).get("conc", [None] * 6)[q] for r in rows]
            cell["conc"][lab] = mean_ci([float(x) if x is not None else float("nan") for x in v],
                                        folds, f"{nm}.conc.{lab}", n_boot=1200)
        cell["conc_near"] = {}
        for q, lab in enumerate(DELTA_LABELS):
            v = [(r["ch"].get(nm, {}).get("conc_near") or [None] * 6)[q] for r in rows]
            cell["conc_near"][lab] = mean_ci([float(x) if x is not None else float("nan")
                                              for x in v], folds, f"{nm}.concnear.{lab}",
                                             n_boot=1200)

        def _res(d):
            for lab in DELTA_LABELS:                  # the SMALLEST bin that stays above 0.5
                c = d[lab]
                if np.isfinite(c.get("mean", float("nan"))) and c["mean"] > 0.5 \
                        and c["fold_ci"][0] > 0.5:
                    return lab
            return None
        cell["resolution"] = _res(cell["conc"])
        cell["resolution_near"] = _res(cell["conc_near"])
        # ---- F-R1, evaluated exactly as registered
        c1a = cell["rho_A_part"]["mean"] >= BAR_RHO and cell["rho_A_part"]["excludes_zero"]
        ctr = cell["contrast_A_minus_ANCHOR"]
        c1b = np.isfinite(ctr.get("effect", np.nan)) and ctr["effect"] >= BAR_RHO_MARGIN and ctr["excludes_zero"]
        c2a = cell["pref_near"]["mean"] >= BAR_PREF and cell["pref_near"]["fold_ci"][0] > 0.5
        cpr = cell["contrast_pref"]
        c2b = np.isfinite(cpr.get("effect", np.nan)) and cpr["effect"] >= BAR_PREF_MARGIN and cpr["excludes_zero"]
        cell["F_R1"] = dict(clause_i=bool(c1a and c1b), clause_i_ordering=bool(c1a),
                            clause_i_margin=bool(c1b), clause_ii=bool(c2a and c2b),
                            clause_ii_pref=bool(c2a), clause_ii_margin=bool(c2b),
                            FIRES=bool(c1a and c1b and c2a and c2b))
        tab[nm] = cell
    o["channels_table"] = tab

    # ---------------- multiplicity: one max-over-channels null per column
    o["A4_multiplicity"] = {}
    for col in ("rho_A_part", "pref_near"):
        M = np.array([[r["ch"].get(nm, {}).get(col, np.nan) for nm in names] for r in rows], float)
        if col.startswith("pref"):
            M = M - 0.5
        o["A4_multiplicity"][col] = signflip_pmax(M)

    # ---------------- the leave-fold-out combination
    core = [nm for nm in names if all(nm in r["ch"] and r["ch"][nm].get("d_near") is not None
                                      for r in rows)]
    o["LFO_combination"] = dict(channels_used=core, **{k: v for k, v in
                                (lfo_combination(rows, core) or {}).items()})
    for tag in ("real", "shuffled"):
        d = o["LFO_combination"].get(tag)
        if isinstance(d, dict) and "_pn" in d:
            d.pop("_pn", None); d.pop("_pc", None); d.pop("folds", None)
    o["LFO_garbage_check"] = dict(**{k: v for k, v in
                                     (lfo_combination(rows, core, key="d_far") or {}).items()})
    for tag in ("real", "shuffled"):
        d = o["LFO_garbage_check"].get(tag)
        if isinstance(d, dict) and "_pn" in d:
            d.pop("_pn", None); d.pop("_pc", None); d.pop("folds", None)

    # ---------------- the verdict
    fired = [nm for nm in names if tab[nm]["F_R1"]["FIRES"]]
    o["VERDICT"] = dict(
        F_R1_FIRES=bool(fired), channels_firing=fired,
        clause_i_only=[nm for nm in names if tab[nm]["F_R1"]["clause_i"] and not tab[nm]["F_R1"]["clause_ii"]],
        clause_ii_only=[nm for nm in names if tab[nm]["F_R1"]["clause_ii"] and not tab[nm]["F_R1"]["clause_i"]],
        best_rho_A_part=max(names, key=lambda nm: tab[nm]["rho_A_part"]["mean"]
                            if np.isfinite(tab[nm]["rho_A_part"]["mean"]) else -9),
        best_pref_near=max(names, key=lambda nm: tab[nm]["pref_near"]["mean"]
                           if np.isfinite(tab[nm]["pref_near"]["mean"]) else -9))

    os.makedirs(RESULTS, exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(o, fh, indent=1, default=float)
    render(o)
    print(f"\nwritten -> {a.out}")
    return o


def render(o):
    t = o["channels_table"]
    ins = o["instrument"]
    print("\n=== INSTRUMENT (every ladder sentence carries this) ===")
    print(f"  torsion-rebuild FLOOR      {ins['rebuild_floor']['mean']:.3f} A  (not 0 A)")
    print(f"  near-native band mean      {ins['near_rmsd']:.3f} A over {ins['n_near']:.1f} rungs/target"
          f"   fallbacks {ins['n_fallback']}")
    print(f"  PROD chain / circ_best     {ins['prod_rmsd']:.3f} / {ins['circbest_rmsd']:.3f} A    "
          f"pool best {ins['pool_best']:.3f} A")
    print(f"  ladder share of the set    {ins['contamination']:.3f}  "
          f"(set-referenced channels use the POOL ONLY: {ins['pool_ref_channels']})")
    print("\n=== A2 realism flatness of ladder A within m (ORACLE; a positive = residual confound) ===")
    for k, v in o["A2_realism_flatness"].items():
        print(f"  {k:8s} {v['mean']:+.3f}  fold CI [{v['fold_ci'][0]:+.3f}, {v['fold_ci'][1]:+.3f}]")
    print("\n=== SI audit: the size-matched twin must zero the PURE functions of Rg ===")
    for k, v in o["SI_audit"].items():
        print(f"  {k:8s} rel sd  raw {v['rel_sd_raw']:.4f} -> SI {v['rel_sd_SI']:.2e}")
    d3 = o["D3_anchor_confound"]
    print(f"\n=== D3 anchor confound  corr(d-to-anchor, d-to-native) within m ===")
    print(f"  median {d3['median']:+.3f}  mean {d3['mean']:+.3f}  [p10 {d3['p10']:+.3f}, p90 {d3['p90']:+.3f}]"
          f"  share > 0.5 {d3['share_above_0p5']:.2f}   CAVEAT TRIGGERS: {d3['caveat_triggers']}")
    d1 = o["D1_locality"]
    print("\n=== D1 locality decomposition (held-out R2 of RMSD-to-native on ladder A) ===")
    for k in ("r2_m_only", "r2_local", "r2_global", "d_local", "d_global"):
        v = d1[k]
        print(f"  {k:11s} {v['mean']:+.4f}  fold CI [{v['fold_ci'][0]:+.4f}, {v['fold_ci'][1]:+.4f}]")

    print("\n=== THE TABLE (every channel; nothing selected) ===")
    hdr = (f"{'channel':16s} {'rhoA_part':>10s} {'rhoANCH_p':>10s} {'A-ANCH':>8s} {'B-ANCH':>8s} "
           f"{'rhoB_pa':>8s} {'rho_rg':>7s} {'prefN':>6s} {'prefPool':>8s} {'pref-ctl':>8s} "
           f"{'pctNat':>7s} {'res_all':>9s} {'res_near':>9s}  F-R1")
    print(hdr); print("-" * len(hdr))
    order = sorted(t, key=lambda nm: -(t[nm]["rho_A_part"]["mean"]
                                       if np.isfinite(t[nm]["rho_A_part"]["mean"]) else -9))
    for nm in order:
        c = t[nm]
        f = c["F_R1"]
        flag = "FIRES" if f["FIRES"] else ("i" if f["clause_i"] else "") + ("ii" if f["clause_ii"] else "") or "-"
        print(f"{nm:16s} {c['rho_A_part']['mean']:+10.3f} {c['rho_ANCHOR_part']['mean']:+10.3f} "
              f"{c['contrast_A_minus_ANCHOR'].get('effect', float('nan')):+8.3f} "
              f"{c['contrast_B_minus_ANCHOR'].get('effect', float('nan')):+8.3f} "
              f"{c['rho_B_part_anchor']['mean']:+8.3f} {c['rho_rg_A']['mean']:+7.3f} "
              f"{c['pref_near']['mean']:6.3f} {c['pref_pool']['mean']:8.3f} "
              f"{c['contrast_pref'].get('effect', float('nan')):+8.3f} "
              f"{c['pctile_nat_in_A']['mean']:7.3f} {str(c['resolution']):>9s} "
              f"{str(c['resolution_near']):>9s}  {flag}")

    print("\n=== A4 multiplicity (max over channels, per-target sign-flip null, 500 draws) ===")
    for k, v in o["A4_multiplicity"].items():
        print(f"  {k:12s} observed max {v['observed_max']:.3f}  null mean {v['null_mean']:.3f} "
              f"p95 {v['null_p95']:.3f}  p_max {v['p_max']:.3f}")

    print("\n=== LFO combination (fit on 4 folds, evaluated on the 5th) ===")
    for tag in ("real", "shuffled"):
        d = o["LFO_combination"].get(tag)
        if isinstance(d, dict) and "pref_near" in d:
            print(f"  near-vs-PROD  {tag:9s} pref_near {d['pref_near']:.3f}  "
                  f"pool-member control {d['pref_pool']:.3f}  margin {d['margin']:+.3f}")
    for tag in ("real",):
        d = o["LFO_garbage_check"].get(tag)
        if isinstance(d, dict) and "pref_near" in d:
            print(f"  GARBAGE CHECK far(>=3A)-vs-PROD    pref {d['pref_near']:.3f}  "
                  f"control {d['pref_pool']:.3f}  margin {d['margin']:+.3f}")

    v = o["VERDICT"]
    print("\n=== VERDICT ON F-R1 ===")
    print(f"  FIRES: {v['F_R1_FIRES']}   channels firing: {v['channels_firing'] or 'NONE'}")
    print(f"  clause (i) only (ordering, no preference): {v['clause_i_only'] or 'none'}")
    print(f"  clause (ii) only (preference, no ordering): {v['clause_ii_only'] or 'none'}")
    print(f"  best partialled rho_A: {v['best_rho_A_part']}   best pref_near: {v['best_pref_near']}")


if __name__ == "__main__":
    main()
