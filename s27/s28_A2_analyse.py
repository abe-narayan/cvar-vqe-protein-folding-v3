#!/usr/bin/env python
"""s27/s28_A2_analyse.py -- statistics for S28 lane A2 (`s27/PREREG_S28_A.md` ADDENDUM 4).

Reads `s27/results/s28_A2_cosine_rows.jsonl` (ORACLE diagnostic), `s28_A2_ladder_rows.jsonl`
(point cloud) and the A2 arms inside `s28_A_chain_rows.jsonl` (built chain); prints every
contrast with `ST.fmt`; writes `s27/results/s28_A2_summary.json`.
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

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A_analyse as AN        # noqa: E402
from s27 import s28_A2_local as L          # noqa: E402

RESULTS = AN.RESULTS
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def load(name):
    f = os.path.join(RESULTS, name)
    rows = {}
    if os.path.exists(f):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line); rows[r["pdb"]] = r
                except Exception:
                    pass
    return rows


def sign_summary(x, label):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    n = len(x)
    if n == 0:
        say("  %-40s n=0" % label); return None
    se = x.std(ddof=1) / np.sqrt(n) if n > 1 else float("nan")
    pos = int((x > 0).sum())
    from math import comb
    # two-sided exact sign test against P(positive) = 1/2
    p = min(1.0, 2 * sum(comb(n, k) for k in range(min(pos, n - pos) + 1)) / 2 ** n)
    say("  %-40s mean %+.3f  SE %.3f  median %+.3f  positive %d/%d  sign-test p %.3g" % (label, x.mean(), se, np.median(x), pos, n, p))
    return dict(mean=float(x.mean()), se=float(se), median=float(np.median(x)), n_pos=pos, n=n, p_sign=float(p))


def main():
    pdbs = sorted(t["pdb"] for t in I.targets())
    folds = ST.pinned_folds(pdbs)
    fail = np.array([p in I.FAIL18 for p in pdbs])
    cos, lad, chn = load("s28_A2_cosine_rows.jsonl"), load("s28_A2_ladder_rows.jsonl"), AN.load_rows("chain")
    summary = dict(n_cosine=len(cos), n_ladder=len(lad), n_chain=len(chn), contrasts={})
    say("S28 lane A2 analysis: cosine rows %d, ladder rows %d, chain rows %d" % (len(cos), len(lad), len(chn)))

    if cos:
        say("\n" + "=" * 100 + "\nA2.1 ORACLE DIAGNOSTIC: cos(-grad, u) at the production point, rigid-body removed  (u is the ORACLE direction to the native)  n=%d" % len(cos))
        chans = list(next(iter(cos.values()))["cos"].keys())
        summary["cosine"] = {}
        for ch in chans:
            v = np.array([cos[p]["cos"][ch] if p in cos else np.nan for p in pdbs])
            lab = ch + (" (smoothed FD, step function)" if ch in L.CHANNELS_FD else " (analytic)")
            say("\n  %s" % lab)
            summary["cosine"][ch] = dict(all=sign_summary(v, "all 126"),
                                         fail18=sign_summary(v[fail], "FAIL18"),
                                         other=sign_summary(v[~fail], "other 108"))
        ref = np.array([np.mean(np.abs(cos[p]["cos_random_ref"])) if p in cos else np.nan for p in pdbs])
        say("\n  random-direction reference: mean |cos| of a random shape field with u = %.3f (the scale of a meaningless cosine at 3n-6 dof)" % np.nanmean(ref))
        summary["cos_random_ref_mean_abs"] = float(np.nanmean(ref))
        gr = np.array([cos[p]["grad_rms"]["DIS"] if p in cos else np.nan for p in pdbs])
        ur = np.array([cos[p]["u_rms"] if p in cos else np.nan for p in pdbs])
        say("  |grad S~| RMS per atom: mean %.4f A^-1; |u| RMS (distance to the native in the frame): mean %.3f A" % (np.nanmean(gr), np.nanmean(ur)))
        # correlation of the DIS cosine with the production RMSD
        from scipy.stats import spearmanr
        v = np.array([cos[p]["cos"]["DIS"] if p in cos else np.nan for p in pdbs]); r0 = np.array([cos[p]["rmsd_prod_cloud"] if p in cos else np.nan for p in pdbs])
        ok = np.isfinite(v) & np.isfinite(r0)
        say("  Spearman(cos_DIS, production RMSD) = %.3f (n=%d)" % (spearmanr(v[ok], r0[ok]).correlation, ok.sum()))

    if lad:
        say("\n" + "=" * 100 + "\nA2.2 DEPLOYABLE STEP LADDER, POINT CLOUD (intermediate; vs production; negative = better)  n=%d" % len(lad))
        pc = np.array([lad[p]["arms"]["prod"]["rmsd_cloud"] if p in lad else np.nan for p in pdbs])
        summary["ladder_cloud"] = {}
        for e in L.STEPS:
            say("")
            v = np.array([lad[p]["arms"]["step_e%g" % e]["rmsd_cloud"] if p in lad else np.nan for p in pdbs])
            r = AN.contrast(v, pc, folds, pdbs, "gradient step e=%g A vs production (point cloud)" % e)
            summary["contrasts"]["cloud:step_e%g" % e] = r
            R = np.column_stack([np.array([lad[p]["arms"]["rand%d_e%g" % (d, e)]["rmsd_cloud"] if p in lad else np.nan for p in pdbs]) for d in range(L.N_RAND)])
            rm = np.nanmean(R, 1)
            say("")
            summary["contrasts"]["cloud:rand_mean8_e%g" % e] = AN.contrast(rm, pc, folds, pdbs, "random direction e=%g A, MEAN of 8, vs production (point cloud)" % e)
            say("")
            summary["contrasts"]["cloud:step_vs_rand_e%g" % e] = AN.contrast(v, rm, folds, pdbs, "gradient step e=%g vs random direction mean-of-8 (point cloud)" % e)
            say("    random best-of-8 (order statistic) mean %.4f" % np.nanmean(np.nanmin(R, 1)))
            cv = np.array([lad[p]["arms"]["circ_e%g" % e]["rmsd_cloud"] if p in lad else np.nan for p in pdbs])
            cP = np.array([lad[p]["arms"]["circP"]["rmsd_cloud"] if p in lad else np.nan for p in pdbs])
            say("")
            summary["contrasts"]["cloud:circ_e%g_vs_circP" % e] = AN.contrast(cv, cP, folds, pdbs, "circuit one-step e=%g vs the family's nearest point to production (point cloud)" % e)
            say("")
            summary["contrasts"]["cloud:circ_e%g_vs_prod" % e] = AN.contrast(cv, pc, folds, pdbs, "circuit one-step e=%g vs production (point cloud)" % e)
            if fail.any():
                d = v - pc
                say("    FAIL18 mean d %+.4f (n=%d)   other 108 mean d %+.4f (n=%d)" % (np.nanmean(d[fail]), fail.sum(), np.nanmean(d[~fail]), (~fail).sum()))
            for nm, vec in (("step", v), ("rand_mean8", rm), ("circ", cv)):
                ss = np.array([lad[p]["arms"][("step_e%g" if nm == "step" else "circ_e%g") % e]["S_smooth"] if p in lad and nm != "rand_mean8" else np.nan for p in pdbs])
                summary["ladder_cloud"]["%s_e%g" % (nm, e)] = float(np.nanmean(vec))
            sp = np.array([lad[p]["S_prod"] if p in lad else np.nan for p in pdbs])
            ss = np.array([lad[p]["arms"]["step_e%g" % e]["S_smooth"] if p in lad else np.nan for p in pdbs])
            say("    S~ production %.4f -> after the step %.4f (mean); the step lowers S~ on %d/%d" % (np.nanmean(sp), np.nanmean(ss), int(np.nansum(ss < sp)), int(np.isfinite(ss).sum())))
        say("")
        cP = np.array([lad[p]["arms"]["circP"]["rmsd_cloud"] if p in lad else np.nan for p in pdbs])
        res = np.array([lad[p]["arms"]["circP"]["residual_rms_to_prod"] if p in lad else np.nan for p in pdbs])
        summary["contrasts"]["cloud:circP_vs_prod"] = AN.contrast(cP, pc, folds, pdbs, "the circuit family's nearest point to production vs production (point cloud)")
        say("    residual RMS of the family's nearest point to C0: mean %.3f A, median %.3f, max %.3f" % (np.nanmean(res), np.nanmedian(res), np.nanmax(res)))
        summary["circP_residual_rms"] = dict(mean=float(np.nanmean(res)), median=float(np.nanmedian(res)), max=float(np.nanmax(res)))

    a2_arms = sorted({a for r in chn.values() for a in r["arms"] if a.startswith(("step_", "rand", "circ_e", "circP"))})
    if chn and a2_arms:
        say("\n" + "=" * 100 + "\nA2.2 BUILT CHAIN (the verdict basis; vs production re-projected in the same job)  arms: %s" % a2_arms)
        pch = AN.arm_vec(chn, pdbs, "prod", "rmsd_chain")
        for arm in a2_arms:
            v = AN.arm_vec(chn, pdbs, arm, "rmsd_chain")
            say("")
            r = AN.contrast(v, pch, folds, pdbs, "%s vs production (BUILT CHAIN)" % arm)
            summary["contrasts"]["chain:" + arm] = r
            if r is not None:
                ok = np.isfinite(v) & np.isfinite(pch); d = v - pch
                say("    FAIL18 mean d %+.4f (n=%d)   other 108 mean d %+.4f (n=%d)" % (np.nanmean(d[fail & ok]), (fail & ok).sum(), np.nanmean(d[~fail & ok]), (~fail & ok).sum()))
        for e in L.STEPS:
            s_arm, r_arms = "step_e%g" % e, [a for a in a2_arms if a.startswith("rand") and a.endswith("_e%g" % e)]
            if s_arm in a2_arms and r_arms:
                rm = np.nanmean(np.column_stack([AN.arm_vec(chn, pdbs, a, "rmsd_chain") for a in r_arms]), 1)
                say("")
                summary["contrasts"]["chain:step_vs_rand_e%g" % e] = AN.contrast(AN.arm_vec(chn, pdbs, s_arm, "rmsd_chain"), rm, folds, pdbs, "gradient step e=%g vs random direction (mean of %d draws) (BUILT CHAIN)" % (e, len(r_arms)))

    path = os.path.join(RESULTS, "s28_A2_summary.json")
    ST.save_atomic(path, dict(summary, text="\n".join(OUT)), module_file=__file__)
    say("\nwritten:", path)


if __name__ == "__main__":
    main()
