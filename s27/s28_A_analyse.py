#!/usr/bin/env python
"""s27/s28_A_analyse.py -- statistics for S28 lane A (`s27/PREREG_S28_A.md` section 7).

Reads `s27/results/s28_A_{oracle,recog,chain}_rows.jsonl`, prints every contrast with
`s24.stats_lib.fmt` verbatim, and writes `s27/results/s28_A_summary.json` (provenance-stamped).
Every contrast is paired per target; d = arm - production; NEGATIVE is better.  Bases are never
mixed: the built-chain table and the point-cloud table are separate blocks.

    python s27/s28_A_analyse.py [--seed 0]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def load_rows(mode, seed=0):
    f = os.path.join(RESULTS, f"s28_A_{mode}_rows.jsonl" if seed == 0 else f"s28_A_{mode}_rows_seed{seed}.jsonl")
    rows = {}
    if os.path.exists(f):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line); rows[r["pdb"]] = r
                except Exception:
                    pass
    return rows


def arm_vec(rows, pdbs, arm, key):
    return np.array([rows[p]["arms"].get(arm, {}).get(key, np.nan) if p in rows else np.nan for p in pdbs], float)


def zero_info_vec(rows, pdbs, key):
    """PREREG addendum 3 (e): the target's zero-information value, the mean over the untrained draws
    (point cloud) or draw 0 (chain), used to score an UNDEFINED readout; never a repair."""
    out = np.full(len(pdbs), np.nan)
    for i, p in enumerate(pdbs):
        if p not in rows:
            continue
        vals = [v.get(key, np.nan) for a, v in rows[p]["arms"].items() if a.startswith("untr_")]
        vals = [v for v in vals if np.isfinite(v)]
        out[i] = float(np.mean(vals)) if vals else np.nan
    return out


def arm_vec_e(rows, pdbs, arm, key, zero_info):
    """`arm_vec` with rule (e) applied: undefined (NaN) deployable readouts scored at the
    zero-information value.  Returns (vector, number substituted)."""
    v = arm_vec(rows, pdbs, arm, key)
    und = np.array([bool(rows[p]["arms"].get(arm, {}).get("undefined", False)) if p in rows else False for p in pdbs])
    sub = und & np.isfinite(zero_info)
    v = v.copy(); v[sub] = zero_info[sub]
    return v, int(sub.sum())


def contrast(a, b, folds, pdbs, label):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        say("  %s: NOT MEASURED (%d finite pairs)" % (label, ok.sum()))
        return None
    r = ST.compare(a[ok], b[ok], folds[ok], names=[p for p, o in zip(pdbs, ok) if o], label=label)
    if ok.sum() < len(a):
        r["n_dropped_nonfinite"] = int((~ok).sum())
        say("  [%d targets dropped: arm undefined (NaN)]" % (~ok).sum())
    say(ST.fmt(r))
    return r


def nested_lfo(M, folds, pdbs, names):
    """Leave-fold-out choice of a grid column: for each held-out fold, the column with the best
    mean on the other four folds is applied to it.  Returns the held-out vector and the choices."""
    out = np.full(M.shape[0], np.nan)
    choice = {}
    for f in sorted(set(folds.tolist())):
        tr = folds != f
        j = int(np.argmin(np.nanmean(M[tr], 0)))
        out[~tr] = M[~tr, j]
        choice[int(f)] = names[j]
    return out, choice


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    seed = a.seed
    pdbs = sorted(t["pdb"] for t in I.targets())
    folds = ST.pinned_folds(pdbs)
    fail = np.array([p in I.FAIL18 for p in pdbs])
    orc, rec, chn = load_rows("oracle", seed), load_rows("recog", seed), load_rows("chain", seed)
    summary = dict(seed=seed, n_oracle=len(orc), n_recog=len(rec), n_chain=len(chn), contrasts={})
    say("S28 lane A analysis: oracle rows %d, recog rows %d, chain rows %d (seed %d)" % (len(orc), len(rec), len(chn), seed))

    # ---------------------------------------------------------------- anchors
    if orc:
        pc = arm_vec(orc, pdbs, "prod", "rmsd_cloud")
        pd_ = arm_vec(orc, pdbs, "prod", "rmsd_cloud_deployed")
        dev = np.array([orc[p].get("prod_frame_max_dev", np.nan) if p in orc else np.nan for p in pdbs])
        say("\nANCHORS (point cloud): production frame readout mean %.6f, deployed readout mean %.6f, "
            "max |dev| over targets %.2e (S25/S27 anchor 3.048338)" % (np.nanmean(pc), np.nanmean(pd_), np.nanmax(dev)))
        summary["anchor_cloud_prod"] = float(np.nanmean(pc)); summary["anchor_cloud_deployed"] = float(np.nanmean(pd_))
        summary["prod_frame_max_dev"] = float(np.nanmax(dev))
    if chn:
        pch = arm_vec(chn, pdbs, "prod", "rmsd_chain")
        say("ANCHOR (built chain): production mean %.6f over %d targets (S27 anchor 3.2126, `chain_rows.jsonl :: DIS`)"
            % (np.nanmean(pch), np.isfinite(pch).sum()))
        summary["anchor_chain_prod"] = float(np.nanmean(pch))
        # cross-check against S27's DIS chain rows
        s27 = {}
        f27 = os.path.join(RESULTS, "chain_rows.jsonl")
        if os.path.exists(f27):
            with open(f27, encoding="utf-8") as fh:
                for line in fh:
                    r = json.loads(line)
                    if r["config"] == "DIS":
                        s27[r["pdb"]] = r["rmsd_chain"]
            v27 = np.array([s27.get(p, np.nan) for p in pdbs])
            ok = np.isfinite(v27) & np.isfinite(pch)
            say("  vs S27 DIS chain rows on %d targets: max |diff| %.2e" % (ok.sum(), np.abs(v27[ok] - pch[ok]).max()))
            summary["chain_vs_s27_max_diff"] = float(np.abs(v27[ok] - pch[ok]).max())

    # ---------------------------------------------------------------- ORACLE expressivity
    if orc:
        say("\n" + "=" * 100 + "\nORACLE EXPRESSIVITY (diagnostic; theta chosen against the native; NEVER a result)  point cloud, n=%d" % len(orc))
        tab = {}
        for arm, key, lab in (("prod", "rmsd_cloud", "production uniform top-75 average"),
                              ("oracle_circ", "rmsd_cloud", "ORACLE circuit family, best of 5 starts"),
                              ("oracle_circ", "rmsd_mean_starts", "ORACLE circuit family, MEAN over 5 starts"),
                              ("oracle_aff500", "rmsd_cloud", "ORACLE affine-500 least squares (fixed frame)"),
                              ("oracle_aff75", "rmsd_cloud", "ORACLE affine-75 least squares (fixed frame)"),
                              ("oracle_sub", "rmsd_cloud", "ORACLE random 27-dim affine subspace, mean of 8")):
            v = arm_vec(orc, pdbs, arm, key)
            tab[arm + ":" + key] = float(np.nanmean(v))
            say("  %-55s mean %.4f  median %.4f  <2A %.3f  n=%d" % (lab, np.nanmean(v), np.nanmedian(v), np.nanmean(v < 2.0), np.isfinite(v).sum()))
        for lab, key in (("ORACLE best member of top-75", "oracle_top75_best"), ("ORACLE best member of pool", "oracle_pool_best")):
            v = np.array([orc[p].get(key, np.nan) if p in orc else np.nan for p in pdbs])
            tab[key] = float(np.nanmean(v))
            say("  %-55s mean %.4f" % (lab, np.nanmean(v)))
        summary["oracle_cloud"] = tab
        # the best-of-5 order statistic
        per = np.array([orc[p]["arms"]["oracle_circ"]["per_start"] if p in orc else [np.nan] * 5 for p in pdbs], float)
        ok = np.isfinite(per).all(1)
        if ok.sum() > 10:
            w = ST.best_of_k_within(per[ok])
            say("  best-of-5 starts priced (best_of_k_within): observed gain %+.4f, valid null %+.4f (%.0f%%), k_eff %.2f, %s"
                % (w["observed_gain"], w["null_across_targets"], 100 * w["share_accounted"], w["k_eff"], w["verdict"]))
            summary["oracle_circ_bok"] = {k: w[k] for k in ("observed_gain", "null_across_targets", "share_accounted", "k_eff", "split_half", "verdict")}
        say("\n  ORACLE contrasts (point cloud; negative = the ORACLE family is nearer the native):")
        oc = arm_vec(orc, pdbs, "oracle_circ", "rmsd_cloud")
        pc = arm_vec(orc, pdbs, "prod", "rmsd_cloud")
        summary["contrasts"]["ORACLE_circ_vs_prod_cloud"] = contrast(oc, pc, folds, pdbs, "ORACLE circuit ceiling vs production (point cloud)")
        summary["contrasts"]["ORACLE_circ_vs_sub_cloud"] = contrast(oc, arm_vec(orc, pdbs, "oracle_sub", "rmsd_cloud"), folds, pdbs, "ORACLE circuit ceiling vs ORACLE random-27-subspace ceiling (point cloud)")
        summary["contrasts"]["ORACLE_circ_vs_aff75_cloud"] = contrast(oc, arm_vec(orc, pdbs, "oracle_aff75", "rmsd_cloud"), folds, pdbs, "ORACLE circuit ceiling vs ORACLE affine-75 ceiling (point cloud)")
        tb = np.array([orc[p].get("oracle_top75_best", np.nan) if p in orc else np.nan for p in pdbs])
        summary["contrasts"]["ORACLE_circ_vs_top75best_cloud"] = contrast(oc, tb, folds, pdbs, "ORACLE circuit ceiling vs ORACLE best single member of top-75 (point cloud)")
        # sign structure of the ORACLE optimum
        for key, lab in (("frac_neg", "fraction of negative weights"), ("neg_mass", "negative mass"), ("ess", "effective members 1/sum w^2"), ("denom", "denominator sum psi"), ("rg", "Rg"), ("bond", "mean virtual bond")):
            v = arm_vec(orc, pdbs, "oracle_circ", key)
            say("  ORACLE circuit optimum: %-32s mean %.4f median %.4f" % (lab, np.nanmean(v), np.nanmedian(v)))
        say("  pool: mean Rg %.3f, mean bond %.3f; production: Rg %.3f, bond %.3f"
            % (np.nanmean([orc[p]["pool_rg_mean"] for p in orc]), np.nanmean([orc[p]["pool_bond_mean"] for p in orc]),
               np.nanmean(arm_vec(orc, pdbs, "prod", "rg")), np.nanmean(arm_vec(orc, pdbs, "prod", "bond"))))
        # does the objective prefer the ORACLE structure?  S~ at the oracle vs at production
        so = arm_vec(orc, pdbs, "oracle_circ", "S_smooth"); sp = arm_vec(orc, pdbs, "prod", "S_smooth")
        ok = np.isfinite(so) & np.isfinite(sp)
        say("  S~ (native-free objective term) at the ORACLE circuit optimum vs at production: mean %.4f vs %.4f; "
            "oracle lower on %d/%d targets" % (so[ok].mean(), sp[ok].mean(), (so[ok] < sp[ok]).sum(), ok.sum()))
        summary["S_at_oracle_vs_prod"] = dict(oracle=float(so[ok].mean()), prod=float(sp[ok].mean()), n_oracle_lower=int((so[ok] < sp[ok]).sum()), n=int(ok.sum()))

    # ---------------------------------------------------------------- BUILT CHAIN (the verdict basis)
    if chn:
        say("\n" + "=" * 100 + "\nBUILT CHAIN (the reporting basis; every arm vs PRODUCTION; negative = better)  n=%d" % len(chn))
        pch = arm_vec(chn, pdbs, "prod", "rmsd_chain")
        arms = sorted({a for r in chn.values() for a in r["arms"]} - {"prod"})
        order = [a for a in ("circ_l0_i80", "circ_l0.3_i80", "circ_l1_i80", "circ_l3_i80") if a in arms]
        order += [a for a in arms if a not in order and not a.startswith("oracle")]
        order += [a for a in arms if a.startswith("oracle")]
        tab = {}
        zi_chain = arm_vec(chn, pdbs, "untr_0", "rmsd_chain")
        for arm in order:
            v = arm_vec(chn, pdbs, arm, "rmsd_chain")
            n_sub = 0
            if not arm.startswith("oracle"):
                v, n_sub = arm_vec_e(chn, pdbs, arm, "rmsd_chain", zi_chain)
            lab = ("ORACLE " if arm.startswith("oracle") else "") + arm + " vs production (BUILT CHAIN)"
            say("")
            if n_sub:
                say("  [rule (e): %d undefined readouts scored at the zero-information value]" % n_sub)
            r = contrast(v, pch, folds, pdbs, lab)
            summary["contrasts"]["chain:" + arm] = r
            tab[arm] = float(np.nanmean(v))
            if r is not None and not arm.startswith("oracle"):
                ok = np.isfinite(v) & np.isfinite(pch)
                d = v - pch
                say("    FAIL18: mean d %+.4f (n=%d)   other 108: mean d %+.4f (n=%d)"
                    % (np.nanmean(d[fail & ok]), (fail & ok).sum(), np.nanmean(d[~fail & ok]), (~fail & ok).sum()))
        summary["chain_means"] = tab
        # nested leave-fold-out choice of lam among the circuit arms
        grid = [a for a in ("circ_l0.3_i80", "circ_l1_i80", "circ_l3_i80") if a in arms]
        if len(grid) == 3:
            M = np.column_stack([arm_vec(chn, pdbs, a, "rmsd_chain") for a in grid])
            ok = np.isfinite(M).all(1) & np.isfinite(pch)
            if ok.sum() > 20:
                held, choice = nested_lfo(M[ok], folds[ok], [p for p, o in zip(pdbs, ok) if o], grid)
                say("\n  lam chosen leave-fold-out (built chain): choices per held-out fold %s" % choice)
                summary["contrasts"]["chain:circ_lam_LFO"] = contrast(held, pch[ok], folds[ok], [p for p, o in zip(pdbs, ok) if o], "circuit, lam chosen leave-fold-out, HELD-OUT vs production (BUILT CHAIN)")
                w = ST.best_of_k_within(M[ok])
                say("  per-target min over the lam grid priced: observed %+.4f, valid null %+.4f (%.0f%%), k_eff %.2f, split-half %+.4f, %s"
                    % (w["observed_gain"], w["null_across_targets"], 100 * w["share_accounted"], w["k_eff"], w["split_half"], w["verdict"]))
                summary["chain_lam_grid_bok"] = {k: w[k] for k in ("observed_gain", "null_across_targets", "share_accounted", "k_eff", "split_half", "verdict")}
        # the circuit against its classical controls on the chain (F2)
        for arm in ("circ_l0.3_i80", "circ_l1_i80", "circ_l3_i80"):
            lam = arm.split("_")[1][1:]
            for ctrl in ("a500_rand_m_l" + lam, "a75_rand_m_l" + lam, "sub0_rand_m_l" + lam, "simplex_rand_m_l" + lam):
                if arm in arms and ctrl in arms:
                    say("")
                    summary["contrasts"]["chain:%s_vs_%s" % (arm, ctrl)] = contrast(arm_vec(chn, pdbs, arm, "rmsd_chain"), arm_vec(chn, pdbs, ctrl, "rmsd_chain"), folds, pdbs, "%s vs %s (BUILT CHAIN)" % (arm, ctrl))

    # ---------------------------------------------------------------- POINT CLOUD (intermediate)
    if rec:
        say("\n" + "=" * 100 + "\nPOINT CLOUD (intermediate; every arm vs production point cloud; negative = better)  n=%d" % len(rec))
        pc = arm_vec(orc, pdbs, "prod", "rmsd_cloud") if orc else np.full(len(pdbs), np.nan)
        arms = sorted({a for r in rec.values() for a in r["arms"]})
        tab = {}
        zi_cloud = zero_info_vec(rec, pdbs, "rmsd_cloud")
        summary["n_undefined"] = {a: int(sum(bool(rec[p]["arms"].get(a, {}).get("undefined", False)) for p in rec)) for a in arms}
        say("  undefined readouts per arm (rule (e)): %s" % {a: n for a, n in summary["n_undefined"].items() if n})
        # circuit arms first
        for arm in [a for a in arms if a.startswith("circ") or a.startswith("diag")]:
            v, n_sub = arm_vec_e(rec, pdbs, arm, "rmsd_cloud", zi_cloud)
            say("")
            if n_sub:
                say("  [rule (e): %d undefined readouts scored at the zero-information value]" % n_sub)
            r = contrast(v, pc, folds, pdbs, arm + " vs production (point cloud)")
            summary["contrasts"]["cloud:" + arm] = r
            tab[arm] = dict(mean=float(np.nanmean(v)), median=float(np.nanmedian(v)), n_undefined=int(np.isnan(v).sum()))
            for key in ("frac_neg", "neg_mass", "ess", "denom", "pad_mass", "rg", "bond", "S_shipped", "S_smooth", "cvar", "H", "F"):
                tab[arm][key] = float(np.nanmean(arm_vec(rec, pdbs, arm, key)))
            say("    diag: frac_neg %.3f  neg_mass %.3f  ess %.1f  |denom| med %.3f  pad %.3f  Rg %.2f  bond %.2f  S %.3f  S~ %.3f  cvar %.3f  H %.3f  F %.3f  undefined %d"
                % (tab[arm]["frac_neg"], tab[arm]["neg_mass"], tab[arm]["ess"], np.nanmedian(np.abs(arm_vec(rec, pdbs, arm, "denom"))), tab[arm]["pad_mass"],
                   tab[arm]["rg"], tab[arm]["bond"], tab[arm]["S_shipped"], tab[arm]["S_smooth"], tab[arm]["cvar"], tab[arm]["H"], tab[arm]["F"], tab[arm]["n_undefined"]))
            if arm == "circ_l0_i80":
                g = arm_vec(rec, pdbs, arm, "gate_lam0_max_abs_dp")
                say("    lam=0 soundness gate: max |p - run_cvar_vqe p| over targets %.2e" % np.nanmax(g))
                summary["gate_lam0_max_abs_dp"] = float(np.nanmax(g))
        # classical controls: a compact table of means, then the contrasts that matter
        say("\n  classical controls, point cloud means (arm: mean / median / n undefined):")
        for arm in [a for a in arms if not (a.startswith("circ") or a.startswith("diag") or a.startswith("untr"))]:
            v = arm_vec(rec, pdbs, arm, "rmsd_cloud")
            tab[arm] = dict(mean=float(np.nanmean(v)), median=float(np.nanmedian(v)), n_undefined=int(np.isnan(v).sum()),
                            frac_neg=float(np.nanmean(arm_vec(rec, pdbs, arm, "frac_neg"))), ess=float(np.nanmean(arm_vec(rec, pdbs, arm, "ess"))),
                            rg=float(np.nanmean(arm_vec(rec, pdbs, arm, "rg"))), S_smooth=float(np.nanmean(arm_vec(rec, pdbs, arm, "S_smooth"))))
            say("    %-28s %.4f / %.4f / %d   frac_neg %.2f  ess %.1f  Rg %.2f  S~ %.3f" % (arm, tab[arm]["mean"], tab[arm]["median"], tab[arm]["n_undefined"], tab[arm]["frac_neg"], tab[arm]["ess"], tab[arm]["rg"], tab[arm]["S_smooth"]))
        # the subspace mean over 8 (the (4b) control) per lam and budget
        for lam in ("0.3", "1", "3"):
            for bud in ("m", "c"):
                cols = [a for a in arms if a.startswith("sub") and a.endswith("_rand_%s_l%s" % (bud, lam))]
                if len(cols) == 8:
                    Ms = np.column_stack([arm_vec(rec, pdbs, a, "rmsd_cloud") for a in cols])
                    tab["sub_mean8_rand_%s_l%s" % (bud, lam)] = float(np.nanmean(np.nanmean(Ms, 1)))
                    tab["sub_best8_rand_%s_l%s" % (bud, lam)] = float(np.nanmean(np.nanmin(Ms, 1)))
                    say("    sub mean-of-8 rand %s lam %s: %.4f   (best-of-8, an order statistic: %.4f)" % (bud, lam, tab["sub_mean8_rand_%s_l%s" % (bud, lam)], tab["sub_best8_rand_%s_l%s" % (bud, lam)]))
        # untrained
        cols = [a for a in arms if a.startswith("untr_")]
        if cols:
            Mu = np.column_stack([arm_vec(rec, pdbs, a, "rmsd_cloud") for a in cols])
            tab["untrained_mean16"] = float(np.nanmean(np.nanmean(Mu, 1))); tab["untrained_best16"] = float(np.nanmean(np.nanmin(Mu, 1)))
            say("    untrained circuit read as signed weights: mean over 16 draws %.4f; best-of-16 (order statistic) %.4f; undefined %d/%d"
                % (tab["untrained_mean16"], tab["untrained_best16"], int(np.isnan(Mu).sum()), Mu.size))
            summary["contrasts"]["cloud:untrained_mean16"] = contrast(np.nanmean(Mu, 1), pc, folds, pdbs, "untrained circuit (mean of 16 draws) vs production (point cloud)")
        summary["cloud"] = tab
        # F2 on the point cloud: circuit vs the matched controls
        for arm in ("circ_l0.3_i80", "circ_l1_i80", "circ_l3_i80"):
            lam = arm.split("_")[1][1:]
            cols = [a for a in arms if a.startswith("sub") and a.endswith("_rand_m_l" + lam)]
            if arm in arms and len(cols) == 8:
                say("")
                summary["contrasts"]["cloud:%s_vs_sub_mean8" % arm] = contrast(arm_vec(rec, pdbs, arm, "rmsd_cloud"), np.nanmean(np.column_stack([arm_vec(rec, pdbs, a, "rmsd_cloud") for a in cols]), 1), folds, pdbs, "%s vs random-27-subspace control (mean of 8, matched budget) (point cloud)" % arm)
            for ctrl in ("a500_rand_m_l" + lam, "a500_rand_c_l" + lam, "simplex_rand_m_l" + lam, "sonly_a500_rand_c"):
                if arm in arms and ctrl in arms:
                    say("")
                    summary["contrasts"]["cloud:%s_vs_%s" % (arm, ctrl)] = contrast(arm_vec(rec, pdbs, arm, "rmsd_cloud"), arm_vec(rec, pdbs, ctrl, "rmsd_cloud"), folds, pdbs, "%s vs %s (point cloud)" % (arm, ctrl))
        # the three-way split, per the coordinator: objective optimum (converged unconstrained) vs what the circuit reached vs emitted
        if "diag_circ_l1_i400" in arms:
            say("")
            summary["contrasts"]["cloud:circ_l1_i80_vs_i400"] = contrast(arm_vec(rec, pdbs, "circ_l1_i80", "rmsd_cloud"), arm_vec(rec, pdbs, "diag_circ_l1_i400", "rmsd_cloud"), folds, pdbs, "circuit lam 1 at 80 iterations vs 400 iterations (convergence diagnostic, point cloud)")
            F80 = arm_vec(rec, pdbs, "circ_l1_i80", "F"); F400 = arm_vec(rec, pdbs, "diag_circ_l1_i400", "F")
            say("    F at 80 vs 400 iterations: mean %.4f vs %.4f; F400 < F80 on %d/%d" % (np.nanmean(F80), np.nanmean(F400), int(np.nansum(F400 < F80)), int(np.isfinite(F80 & 1).sum() if False else np.isfinite(F80).sum())))
        # correlation across targets: ORACLE ceiling vs recognition arm
        if orc and "circ_l1_i80" in arms:
            from scipy.stats import spearmanr
            oc = arm_vec(orc, pdbs, "oracle_circ", "rmsd_cloud"); rc = arm_vec(rec, pdbs, "circ_l1_i80", "rmsd_cloud")
            ok = np.isfinite(oc) & np.isfinite(rc)
            say("\n  ORACLE ceiling vs recognition (circ lam 1) across targets: Spearman %.3f (n=%d)" % (spearmanr(oc[ok], rc[ok]).correlation, ok.sum()))
            so = arm_vec(orc, pdbs, "oracle_circ", "S_smooth"); sa = arm_vec(rec, pdbs, "circ_l1_i80", "S_smooth")
            ok = np.isfinite(so) & np.isfinite(sa)
            say("  S~ at the ORACLE optimum vs at the lam-1 arm's optimum: mean %.4f vs %.4f; ORACLE lower on %d/%d  (if the arm's S~ is lower, the objective is wrong, not the optimiser)"
                % (so[ok].mean(), sa[ok].mean(), (so[ok] < sa[ok]).sum(), ok.sum()))
            summary["S_oracle_vs_arm"] = dict(oracle=float(so[ok].mean()), arm=float(sa[ok].mean()), n_oracle_lower=int((so[ok] < sa[ok]).sum()), n=int(ok.sum()))

    path = os.path.join(RESULTS, "s28_A_summary.json" if seed == 0 else "s28_A_summary_seed%d.json" % seed)
    ST.save_atomic(path, dict(summary, text="\n".join(OUT)), module_file=__file__)
    say("\nwritten:", path)


if __name__ == "__main__":
    main()
