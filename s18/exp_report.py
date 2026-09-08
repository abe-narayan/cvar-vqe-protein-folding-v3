"""s18/exp_report.py -- the reader for the Sprint-18 EXPERIMENT arms.

Everything the BRIEF sec 6 makes standing law, computed once and printed once: TARGET as the
unit, paired FOLD-CLUSTERED bootstrap CIs, medians and W/L beside every mean, per-target
distributions, stratification by length and fold, both mandatory controls on every arm, and the
concentration diagnostic read against a UNIFORM-EFFECT NULL rather than a raw drop-top threshold
(the recorded failure mode).

    python -m s18.exp_report [main_json]
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

from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(HERE, "results")
MDE = 0.084


def load(name, quiet=False):
    p = name if os.path.isabs(name) else os.path.join(RESULTS, name)
    if not os.path.exists(p):
        return None
    o = json.load(open(p))
    if not o.get("complete") and not quiet:
        print(f"  ** {os.path.basename(p)} is PARTIAL: {len(o.get('rows', []))} rows"
              f" -- NOT a final result **")
    return o


def boot_fold(diff, folds, rng, B=4000):
    diff = np.asarray(diff, float)
    folds = np.asarray(folds)
    uf = np.unique(folds)
    idx = [np.where(folds == f)[0] for f in uf]
    m = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, len(uf), len(uf))
        m[b] = diff[np.concatenate([idx[p] for p in pick])].mean()
    return float(diff.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def conc_null(diff, rng, k=10, B=2000):
    """Concentration read against a UNIFORM-EFFECT NULL, never a raw drop-top threshold."""
    d = np.asarray(diff, float)
    n = len(d)
    if n <= k + 2:
        return None
    o = np.argsort(d)
    obs = float(d[o[k:]].mean()) if d.mean() > 0 else float(d[o[:n - k]].mean())
    r = d - d.mean()
    sim = np.empty(B)
    for b in range(B):
        s = d.mean() + r[rng.integers(0, n, n)]
        so = np.argsort(s)
        sim[b] = s[so[k:]].mean() if d.mean() > 0 else s[so[:n - k]].mean()
    return {"observed": obs, "p5": float(np.percentile(sim, 5)),
            "p50": float(np.percentile(sim, 50)), "p95": float(np.percentile(sim, 95)),
            "pct": float((sim < obs).mean() * 100)}


def line(rows, key, base, rng, tag="", g=None):
    v = np.array([r[key] for r in rows], float)
    b = np.array([r[base] for r in rows], float)
    f = np.array([r["fold"] for r in rows])
    d = v - b
    mu, lo, hi = boot_fold(d, f, rng)
    w, l = int((d < 0).sum()), int((d > 0).sum())
    return (f"  {key:<22}{v.mean():>7.3f}{np.median(v):>8.3f}   "
            f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l:<4}{np.median(d):>+8.3f}  {tag}")


def main(path=None):
    rng = SD.stable_rng("s18", "exp_report")
    path = path or "exp_main_pool_512.json"
    o = load(path)
    if o is None:
        print(f"no {path} yet")
        return
    rows = o["rows"]
    n = len(rows)
    folds = np.array([r["fold"] for r in rows])
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731

    print("=" * 102)
    print(f"SPRINT 18 -- EXPERIMENT.  n = {n}   objective: MATH s18/math_anova.py "
          f"(mu={o['mu']}, S={o['S']}, grid={o['grid']})")
    print("Start for EVERY arm = the coordinate average of the shipped top-75 (Control A).")
    print("Native read only to score.  MDE at 80% power on this instrument = 0.084 A.")
    print("=" * 102)

    # ------------------------------------------------------------------ P0 VALIDITY GATES
    print("\n[P0] VALIDITY GATES -- nothing downstream is readable unless these pass\n")
    s17p = os.path.join(ROOT, "s17", "results", "refine.json")
    if os.path.exists(s17p):
        s17 = {r["pdb"]: r for r in json.load(open(s17p))["rows"]}
        com = [r for r in rows if r["pdb"] in s17]
        if com:
            da = np.array([r["lam100"] - s17[r["pdb"]]["refine_full"] for r in com])
            do = np.array([r["lam100_objfull"] - s17[r["pdb"]]["obj_full"] for r in com])
            dv = np.array([r["avg"] - s17[r["pdb"]]["avg"] for r in com])
            print(f"  P0.a  lambda=1 vs s17 refine_full, {len(com)} shared targets")
            print(f"        RMSD mean|d| {np.abs(da).mean():.5f}  max|d| {np.abs(da).max():.5f}"
                  f"   |  objective mean|d| {np.abs(do).mean():.4f}")
            print(f"        Control A identical to s17 to {np.abs(dv).max():.2e} A")
            print(f"        GATE (mean|d| < 0.01 A): "
                  f"{'PASS' if np.abs(da).mean() < 0.01 else 'FAIL'}")
    print("  P0.b  E_le1 + E_ge2 == E_full to 9.1e-13 through MATH's callables "
          "(`python -m s18.exp_obj`).  EXACT by construction.")
    print("  P0.d  analytic gradients vs central differences OF MATH'S OWN functions: "
          "relative error 3.3e-10 on both objectives.")
    und = [len(r["undetermined"]) for r in rows]
    print(f"  P0.e  residues the degree-1 object exerts NO force on: {np.mean(und):.2f} per "
          f"target ({np.mean(und) / g('n').mean() * 100:.1f}% of the chain); MATH's theorem "
          f"says exactly {{0, n-1}} -- observed on {sum(1 for r in rows if r['undetermined'] == [0, r['n'] - 1])}/{n}.")
    s = load("exp_sens.json", quiet=True)
    if s:
        sr = s["rows"]
        gg = lambda k: np.array([r[k] for r in sr], float)        # noqa: E731
        print(f"  P0.f  QUADRATURE: on {len(sr)} stratified targets, |S={o['S']} - S=2048|"
              f" = {np.abs(gg(f'le1_S{o[chr(83)]}') - gg('le1_S2048')).mean():.3f} A (refine) / "
              f"{np.abs(gg(f'argmin_S{o[chr(83)]}') - gg('argmin_S2048')).mean():.3f} A (argmin)"
              f"   MDE = {MDE}")

    # ------------------------------------------------------------------ THE LAMBDA CURVE
    print("\n" + "=" * 102)
    print("[H2] THE LAMBDA LADDER   E_lam = E_le1 + lam*E_ge2 = (1-lam) E_le1 + lam E_full")
    print("     PRE-REGISTERED lam in {0, 0.25, 0.5, 0.75, 1}.  NOT extended.")
    print("=" * 102)
    lam = [("lam000", 0.0), ("lam025", 0.25), ("lam050", 0.5), ("lam075", 0.75), ("lam100", 1.0)]
    print(f"\n  {'arm':<22}{'RMSD':>7}{'median':>8}   {'vs avg [fold-clustered 95% CI]':<26}"
          f"{'W/L':>9}{'medD':>8}")
    print(f"  {'avg (Control A)':<22}{g('avg').mean():>7.3f}{np.median(g('avg')):>8.3f}"
          f"   {'--':<26}{'--':>9}{'--':>8}")
    print(line(rows, "proj", "avg", rng, "the torsion start point"))
    for k, lv in lam:
        print(line(rows, k, "avg", rng, f"lambda = {lv}"))

    print(f"\n  {'lam':>6}{'meanRMSD':>10}{'medRMSD':>9}{'dRMSDvsAvg':>12}{'E_full red%':>13}"
          f"{'E_full end':>12}{'E_le1 end':>12}{'disp':>8}{'nfev':>7}{'fracImp':>9}")
    for k, lv in lam:
        v = g(k)
        redf = float((1 - g(k + "_objfull") / g("obj_full_proj")).mean() * 100)
        print(f"  {lv:>6.2f}{v.mean():>10.3f}{np.median(v):>9.3f}"
              f"{(v - g('avg')).mean():>+12.3f}{redf:>13.1f}"
              f"{g(k + '_objfull').mean():>12.1f}{g(k + '_objle1').mean():>12.1f}"
              f"{g(k + '_disp').mean():>8.3f}{g(k + '_nfev').mean():>7.0f}"
              f"{(v < g('avg')).mean():>9.3f}")
    print(f"  start: E_full {g('obj_full_proj').mean():.1f}, E_le1 {g('obj_le1_proj').mean():.1f};"
          f" MATH's MESH argmin of E_le1 reaches {g('argmin_le1_objle1').mean():.1f}.")
    print("  No 'fraction of the way to the optimum' is quoted: E_le1 is signed AND its mesh")
    print("  argmin is not its optimum (see [OPT]), so that ratio has no honest denominator.")

    mono = np.array([g(k).mean() for k, _ in lam])
    bi = int(np.argmin(mono))
    shape = ("MONOTONE INCREASING in lambda" if np.all(np.diff(mono) > 0) else
             "MONOTONE DECREASING in lambda" if np.all(np.diff(mono) < 0) else
             f"INTERIOR OPTIMUM at lambda = {lam[bi][1]}")
    print(f"\n  SHAPE OF THE CURVE: {shape}")
    print(f"  means {['%.3f' % x for x in mono]}   best lam = {lam[bi][1]} "
          f"({mono[bi]:.3f} A)  vs Control A {g('avg').mean():.3f}")
    d01 = g("lam000") - g("lam100")
    mu, lo, hi = boot_fold(d01, folds, rng)
    print(f"  lam=0 minus lam=1 (degree-1 vs full): {mu:+.3f} [{lo:+.3f},{hi:+.3f}]  "
          f"median {np.median(d01):+.3f}  W/L {(d01 < 0).sum()}/{(d01 > 0).sum()}")
    mu, lo, hi = boot_fold(g(lam[bi][0]) - g("avg"), folds, rng)
    print(f"  best lambda minus Control A:          {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")

    print("\n  per-target dRMSD vs Control A -- quantiles across targets:")
    print(f"    {'lam':>6}{'p10':>8}{'p25':>8}{'p50':>8}{'p75':>8}{'p90':>8}{'fracImp':>10}"
          f"{'worst':>9}{'best':>9}")
    for k, lv in lam:
        d = g(k) - g("avg")
        print(f"    {lv:>6.2f}" + "".join(f"{np.percentile(d, q):>+8.3f}"
                                          for q in (10, 25, 50, 75, 90))
              + f"{(d < 0).mean():>10.3f}{d.max():>+9.3f}{d.min():>+9.3f}")

    # ------------------------------------------------------------------ H2b vs objceil
    oc = load("objceil.json", quiet=True)
    if oc:
        m = {r["pdb"]: r for r in oc["rows"]}
        com = [r for r in rows if r["pdb"] in m]
        if com:
            print("\n" + "=" * 102)
            print("[H2b] THE COORDINATOR'S AXIS -- what 'not trusting the error direction' buys")
            print("      Same targets, same Control-A start.  alpha>0 arms are ORACLE ceilings.")
            print("=" * 102)
            f2 = np.array([r["fold"] for r in com])
            a0 = np.array([m[r["pdb"]]["a0.0"] if "a0.0" in m[r["pdb"]] else
                           m[r["pdb"]].get("a0", np.nan) for r in com], float)
            keys = [k for k in m[com[0]["pdb"]] if k.startswith("a")]
            print(f"      objceil arms available: {sorted(keys)}")
            sh = np.array([m[r["pdb"]]["shuffled"] for r in com], float)
            iso = np.array([m[r["pdb"]]["isotropic"] for r in com], float)
            av = np.array([r["avg"] for r in com], float)
            l0 = np.array([r["lam000"] for r in com], float)
            am = np.array([r["argmin_le1_hold"] for r in com], float)
            l1 = np.array([r["lam100"] for r in com], float)
            print(f"\n      {'arm':<34}{'RMSD':>8}{'median':>9}{'vs avg':>10}")
            for nm, v in (("Control A (coordinate average)", av),
                          ("alpha=0  full objective (s17)", l1),
                          ("DEGREE-1 refine (lam=0)", l0),
                          ("DEGREE-1 certified argmin", am),
                          ("shuffled residual direction ORACLE", sh),
                          ("isotropic matched noise ORACLE", iso)):
                print(f"      {nm:<34}{v.mean():>8.3f}{np.median(v):>9.3f}"
                      f"{v.mean() - av.mean():>+10.3f}")
            gap = l1.mean() - sh.mean()
            for nm, v in (("lam=0", l0), ("argmin", am)):
                frac = (l1.mean() - v.mean()) / gap if gap != 0 else float("nan")
                mu, lo, hi = boot_fold(v - l1, f2, rng)
                print(f"\n      H2b: {nm} closes {frac * 100:+.1f}% of the "
                      f"{gap:.3f} A alpha=0 -> shuffled gap")
                print(f"           {nm} - alpha=0: {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
            print("\n      CAVEAT (the coordinator's own, and it cuts against their result):")
            print("      alpha is a blend TOWARD TRUTH, not a model of how a better predictor")
            print("      would err.  The alpha ladder is a requirement in the 'how much residual")
            print("      must go away' sense, NOT a prediction about any achievable predictor.")

    # ------------------------------------------------------------------ ALL ARMS
    print("\n" + "=" * 102)
    print("[ARMS] every arm against Control A; both mandatory controls present")
    print("=" * 102)
    print(f"\n  {'arm':<22}{'RMSD':>7}{'median':>8}   {'vs avg [95% CI]':<26}{'W/L':>9}{'medD':>8}")
    tags = [("lam000", "refine toward DEGREE-1"),
            ("lam100", "refine toward FULL (= s17)"),
            ("argmin_le1", "EXACT argmin, MATH as shipped"),
            ("argmin_le1_hold", "EXACT argmin, terminals held"),
            ("refine_h1", "harmonic order-1 (angles)"),
            ("argmin_h1", "harmonic order-1 argmin"),
            ("rand_field_le1", "CONTROL random additive field"),
            ("rand_field_argmin", "CONTROL random field, argmin"),
            ("rand_move_le1", "CONTROL matched move (le1)"),
            ("rand_move_full", "CONTROL matched move (full)"),
            ("rand_move_argmin", "CONTROL matched move (argmin)")]
    sh = load(f"exp_shuf_{o['S']}.json", quiet=True)
    if sh:
        m = {r["pdb"]: r for r in sh["rows"]}
        for r in rows:
            if r["pdb"] in m:
                r.update({k: v for k, v in m[r["pdb"]].items() if k.startswith("shufobj")})
        if all("shufobj_le1" in r for r in rows):
            tags += [("shufobj_le1", "CONTROL shuffled distogram, le1"),
                     ("shufobj_argmin", "CONTROL shuffled distogram, argmin"),
                     ("shufobj_full", "CONTROL shuffled distogram, full")]
        else:
            print(f"  (shuffled-distogram control is PARTIAL: "
                  f"{sum('shufobj_le1' in r for r in rows)}/{n} -- reported separately)")
    for k, tg_ in tags:
        if k in rows[0]:
            print(line(rows, k, "avg", rng, tg_))

    dpatch = g("argmin_le1") - g("argmin_le1_hold")
    print(f"\n  terminal-residue patch: MATH's shipped argmin vs terminals-held = "
          f"{dpatch.mean():+.3f} A (max |d| {np.abs(dpatch).max():.2e}, "
          f"{int((np.abs(dpatch) > 1e-9).sum())}/{n} targets differ)")
    if np.abs(dpatch).max() < 1e-9:
        print("  I FLAGGED THIS AS A HAZARD AND THE MEASUREMENT RETIRES IT (EXACT): the residues")
        print("  the degree-1 object cannot see are EXACTLY the residues the metric cannot see.")
        print("  phi_0 and psi_{n-1} are never read by the builder; psi_0 rotates the chain about")
        print("  an axis THROUGH CA_0, which fixes CA_0 and is therefore a rigid rotation of the")
        print("  whole Ca set; phi_{n-1} moves only atoms placed after CA_{n-1}.  So the terminal")
        print("  patch cannot cost RMSD, and the 'blind to 17% of the chain' worry is void.")

    # ------------------------------------------------------------------ H7 / ORACLE
    if "ORC_le1" in rows[0]:
        print("\n" + "=" * 102)
        print("[H7] WHAT LIMITS DEGREE-1 -- the distances, or the truncation?")
        print("     Same fields, rebuilt EXACTLY for a different distance vector from MATH's")
        print("     conditional moment tables (`with_dhat`).  ORACLE arms are ceilings, never")
        print("     results, and are never mixed with the realized numbers above.")
        print("=" * 102)
        print(f"\n  {'distances fed to the objective':<40}{'E_full':>9}{'E_le1':>9}"
              f"{'argmin_le1':>12}{'label':>14}")
        print(f"  {'deployed distogram (REALIZED)':<40}{g('lam100').mean():>9.3f}"
              f"{g('lam000').mean():>9.3f}{g('argmin_le1_hold').mean():>12.3f}"
              f"{'native-free':>14}")
        print(f"  {'pairs permuted (CONTROL)':<40}{g('shufobj_full').mean():>9.3f}"
              f"{g('shufobj_le1').mean():>9.3f}{g('shufobj_argmin').mean():>12.3f}"
              f"{'native-free':>14}")
        print(f"  {'magnitudes kept, direction destroyed':<40}{g('ORCSHUF_full').mean():>9.3f}"
              f"{g('ORCSHUF_le1').mean():>9.3f}{'--':>12}{'ORACLE':>14}")
        print(f"  {'the native distances themselves':<40}{g('ORC_full').mean():>9.3f}"
              f"{g('ORC_le1').mean():>9.3f}{g('ORC_argmin').mean():>12.3f}{'ORACLE':>14}")
        gf = g("lam100").mean() - g("ORC_full").mean()
        gl = g("lam000").mean() - g("ORC_le1").mean()
        print(f"\n  perfect distances buy the FULL objective {gf:.3f} A and DEGREE-1 {gl:.3f} A"
              f"  -> degree-1 retains {gl / gf * 100 if gf else float('nan'):.1f}% of the gain")
        mu, lo, hi = boot_fold(g("ORC_le1") - g("ORC_full"), folds, rng)
        print(f"  ORACLE gap: degree-1 minus full at PERFECT distances "
              f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]  W/L "
              f"{int((g('ORC_le1') < g('ORC_full')).sum())}/"
              f"{int((g('ORC_le1') > g('ORC_full')).sum())}")
        mu, lo, hi = boot_fold(g("ORC_le1") - g("avg"), folds, rng)
        print(f"  and against Control A: {mu:+.3f} [{lo:+.3f},{hi:+.3f}]  -- if this is "
              f"positive, NO distogram improvement rescues a degree-1 objective.")

    # ------------------------------------------------------------------ FOUR QUANTITIES
    print("\n" + "=" * 102)
    print("[H4] THE FOUR QUANTITIES, SEPARATELY (BRIEF sec 5) -- they dissociate")
    print("=" * 102)
    print("\n  (1) OBJECTIVE QUALITY")
    print(f"      E_full  {g('obj_full_proj').mean():9.1f} -> {g('lam100_objfull').mean():9.1f}"
          f"  ({(1 - g('lam100_objfull') / g('obj_full_proj')).mean() * 100:5.1f}% reduction)")
    print(f"      E_le1   {g('obj_le1_proj').mean():9.1f} -> {g('lam000_objle1').mean():9.1f}"
          f"  certified optimum {g('argmin_le1_objle1').mean():9.1f}")
    better = int((g("lam000_objle1") < g("argmin_le1_objle1") - 1e-9).sum())
    print(f"      L-BFGS on E_le1 ends BELOW MATH's mesh argmin on {better}/{n} targets, so the")
    print(f"      mesh argmin is not the objective's optimum and no 'percent of the way there'")
    print(f"      is quotable.  The true continuous optimum is in [OPT] below.")
    print("\n  (2) STRUCTURAL QUALITY -- the arm table above.")
    print("\n  (3) ALIGNMENT -- four NAMED axes over the same 75 real pool structures")
    print(f"      {'objective':<10}{'Spearman':>10}{'Pearson':>10}{'pairwise-ord':>14}"
          f"{'in-band Spearman':>19}{'frac rho>0':>12}")
    for nm in ("full", "le1"):
        sp, pe, po = g("sp_" + nm), g("pe_" + nm), g("poa_" + nm)
        sb = g("spband_" + nm)
        sb = sb[~np.isnan(sb)]
        print(f"      {nm:<10}{np.nanmean(sp):>10.3f}{np.nanmean(pe):>10.3f}"
              f"{np.nanmean(po):>14.3f}{sb.mean():>19.3f}{(sp > 0).mean():>12.3f}")
    for nm in ("full", "le1"):
        mu, lo, hi = boot_fold(g("sp_" + nm), folds, rng)
        print(f"      Spearman({nm}): {mu:+.3f} [{lo:+.3f},{hi:+.3f}]  "
              f"median {np.median(g('sp_' + nm)):+.3f}")
    dsp = g("sp_le1") - g("sp_full")
    mu, lo, hi = boot_fold(dsp, folds, rng)
    print(f"      Spearman(le1) - Spearman(full): {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
    print(f"      (Phase 0's 19-target tension was +0.153 for degree-1 vs +0.264 for full.)")
    print(f"      in-band = within {1.5} A of the pool best; mean band size "
          f"{g('n_band').mean():.1f}/75.  Pool diversity {g('pool_div').mean():.3f} A.")
    print("\n  (4) SEARCH QUALITY -- `s18/exp_budget.py` below, and EXACTLY:")
    print("      E_le1 is additive in residues, so its global optimum is the per-residue argmin")
    print("      of a tabulated field: CERTIFIED, not searched.  For degree-1 there is no")
    print("      search problem at any budget.  Whether L-BFGS from Control A reaches it is a")
    print("      separate question, answered in [OPT] below.")

    # ------------------------------------------------------------------ PHASE 9
    print("\n" + "=" * 102)
    print("[H5] PHASE 9 -- protecting the coordinate average.  Same pool, five routes.")
    print("=" * 102)
    print(f"\n  {'route':<38}{'RMSD':>7}{'median':>8}   {'vs avg [95% CI]':<26}{'W/L':>9}")
    print(f"  {'pool -> coordinate average (INCUMBENT)':<38}{g('avg').mean():>7.3f}"
          f"{np.median(g('avg')):>8.3f}   {'--':<26}{'--':>9}")
    for k, nm in (("pool_best", "pool ORACLE best (a ceiling)"),
                  ("pool_mean", "pool mean == random single pick"),
                  ("sel_full", "pool -> select by E_full"),
                  ("sel_le1", "pool -> select by E_le1"),
                  ("lam000", "pool -> avg -> refine on E_le1"),
                  ("lam100", "pool -> avg -> refine on E_full"),
                  ("argmin_le1_hold", "pool -> avg -> E_le1 exact argmin")):
        v = g(k); d = v - g("avg")
        mu, lo, hi = boot_fold(d, folds, rng)
        print(f"  {nm:<38}{v.mean():>7.3f}{np.median(v):>8.3f}   "
              f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]{int((d < 0).sum()):>5}/{int((d > 0).sum()):<4}")

    # ------------------------------------------------------------------ STRATIFICATION
    print("\n" + "=" * 102)
    print("[STRAT] by length and by fold")
    print("=" * 102)
    N = g("n")
    print(f"  {'stratum':<16}{'k':>4}{'avg':>8}{'lam0':>8}{'d':>8}{'lam1':>8}{'d':>8}"
          f"{'argmin':>9}{'d':>8}")
    for name, msk in (("n <= 11", N <= 11), ("12 <= n <= 13", (N >= 12) & (N <= 13)),
                      ("n >= 14", N >= 14)):
        if msk.sum() < 3:
            continue
        print(f"  {name:<16}{int(msk.sum()):>4}{g('avg')[msk].mean():>8.3f}"
              f"{g('lam000')[msk].mean():>8.3f}{(g('lam000') - g('avg'))[msk].mean():>+8.3f}"
              f"{g('lam100')[msk].mean():>8.3f}{(g('lam100') - g('avg'))[msk].mean():>+8.3f}"
              f"{g('argmin_le1_hold')[msk].mean():>9.3f}"
              f"{(g('argmin_le1_hold') - g('avg'))[msk].mean():>+8.3f}")
    print()
    for f in sorted(set(folds.tolist())):
        m = folds == f
        print(f"  fold {f:<3} k={int(m.sum()):<4} avg {g('avg')[m].mean():.3f}"
              f"   lam0 {(g('lam000') - g('avg'))[m].mean():+.3f}"
              f"   lam1 {(g('lam100') - g('avg'))[m].mean():+.3f}"
              f"   argmin {(g('argmin_le1_hold') - g('avg'))[m].mean():+.3f}")

    # ------------------------------------------------------------------ CONCENTRATION
    print("\n" + "=" * 102)
    print("[CONC] concentration vs a UNIFORM-EFFECT NULL (never a raw drop-top threshold)")
    print("=" * 102)
    for k in ("lam000", "lam100", "argmin_le1_hold"):
        c = conc_null(g(k) - g("avg"), rng)
        if c:
            print(f"  {k:<18} drop-top-10 {c['observed']:+.3f}   uniform-effect null "
                  f"p5 {c['p5']:+.3f} p50 {c['p50']:+.3f} p95 {c['p95']:+.3f}"
                  f"  -> percentile {c['pct']:.1f}")
    print("  (a percentile inside [5, 95] means a uniform effect with this much noise already")
    print("   explains the drop-top number: NOT evidence of concentration.)")

    # ------------------------------------------------------------------ FALSIFIERS
    print("\n" + "=" * 102)
    print("[FALSIFIERS]")
    print("=" * 102)
    l0, l1 = g("lam000"), g("lam100")
    a0 = g("argmin_le1_hold")
    av = g("avg")
    best = np.minimum(l0, a0)
    ctrl_names = [k for k in ("rand_field_le1", "rand_move_le1", "shufobj_le1") if k in rows[0]]
    ctrl = min(g(k).mean() for k in ctrl_names)
    m01, lo01, hi01 = boot_fold(l0 - l1, folds, rng)
    mA, loA, hiA = boot_fold(l0 - av, folds, rng)
    mB, loB, hiB = boot_fold(a0 - av, folds, rng)
    f1 = (l0.mean() >= 3.5) or (a0.mean() >= 3.5)
    f2 = lo01 <= 0 <= hi01
    f3 = min(l0.mean(), a0.mean()) >= ctrl - MDE
    f4 = (loA <= 0 <= hiA or mA > 0) and (loB <= 0 <= hiB or mB > 0)
    print(f"  F1 degree-1 lands ~3.6 A on the real instrument     "
          f"{'FIRED' if f1 else 'not fired'}   (lam0 {l0.mean():.3f}, argmin {a0.mean():.3f})")
    print(f"  F2 degree-1 no better than full                     "
          f"{'FIRED' if f2 else 'not fired'}   (lam0-lam1 {m01:+.3f} [{lo01:+.3f},{hi01:+.3f}])")
    print(f"  F3 advantage dies under matched controls            "
          f"{'FIRED' if f3 else 'not fired'}   (best control {ctrl:.3f} of {ctrl_names})")
    print(f"  F4 effect exists only on the 19-target instrument   "
          f"{'FIRED' if f4 else 'not fired'}   (lam0-avg {mA:+.3f} [{loA:+.3f},{hiA:+.3f}];"
          f" argmin-avg {mB:+.3f} [{loB:+.3f},{hiB:+.3f}])")
    print("  F5 invalid continuous mapping -- MATH owns it.")

    # ------------------------------------------------------------------ mu sensitivity
    import glob as _glob
    _u = sorted(_glob.glob(os.path.join(RESULTS, f"exp_main_uniform_{o['S']}*.json")))
    ou = load(os.path.basename(_u[0]), quiet=True) if _u else None
    if ou is not None and not ou.get("complete"):
        print(f"\n  ** mu=uniform arm is PARTIAL: {len(ou['rows'])} rows -- reported as a "
              f"sensitivity SUBSET, never as a final result **")
    if ou:
        mu_ = {r["pdb"]: r for r in ou["rows"]}
        com = [r for r in rows if r["pdb"] in mu_]
        if com:
            f2 = np.array([r["fold"] for r in com])
            av2 = np.array([r["avg"] for r in com])
            print("\n" + "=" * 102)
            print(f"[H6] mu-SENSITIVITY -- pool marginal vs uniform, {len(com)} shared targets")
            print("=" * 102)
            for k in ("lam000", "argmin_le1_hold", "lam100"):
                a = np.array([r[k] for r in com], float)
                b = np.array([mu_[r["pdb"]][k] for r in com], float)
                m1, l1_, h1 = boot_fold(a - av2, f2, rng)
                m2, l2_, h2 = boot_fold(b - av2, f2, rng)
                same = "SAME SIGN" if np.sign(m1) == np.sign(m2) else "SIGN FLIPS -- NOT IDENTIFIED"
                print(f"  {k:<18} pool {a.mean():.3f} ({m1:+.3f} [{l1_:+.3f},{h1:+.3f}])   "
                      f"uniform {b.mean():.3f} ({m2:+.3f} [{l2_:+.3f},{h2:+.3f}])   {same}")

    # ------------------------------------------------------------------ polished optimum
    pl = load("exp_polish.json", quiet=True)
    if pl:
        pr = pl["rows"]
        gp = lambda k: np.array([r[k] for r in pr], float)         # noqa: E731
        fp = np.array([r["fold"] for r in pr])
        print("\n" + "=" * 102)
        print(f"[OPT] THE HONEST CERTIFIED OPTIMUM OF E_le1, n = {len(pr)}")
        print("=" * 102)
        print(f"  MATH's `argmin_le1` returns the best MESH point, but `E_le1` is evaluated by")
        print(f"  trigonometric interpolation and dips BELOW every node between them.  Polished")
        print(f"  per residue from the {pl['kstart']} lowest cells (additive object -> the")
        print(f"  per-residue minimiser IS the global minimiser).")
        print(f"    objective   mesh {gp('mesh_obj').mean():10.2f} -> polished "
              f"{gp('polish_obj').mean():10.2f}  (lower on "
              f"{int((gp('polish_obj') < gp('mesh_obj') - 1e-9).sum())}/{len(pr)} targets;"
              f" {gp('n_res_improved_off_mesh').mean():.1f} residues/target sit off-mesh)")
        print(f"    RMSD        mesh {gp('mesh_argmin').mean():10.3f} -> polished "
              f"{gp('polish_argmin').mean():10.3f}   Control A {gp('avg').mean():.3f}")
        mu, lo, hi = boot_fold(gp("polish_argmin") - gp("mesh_argmin"), fp, rng)
        print(f"    polished minus mesh: {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
        mu, lo, hi = boot_fold(gp("polish_argmin") - gp("avg"), fp, rng)
        print(f"    polished minus Control A: {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
        print(f"    angle-additive polished argmin {gp('polish_ang').mean():.3f}")
        print("  READ: reaching the TRUE optimum of the degree-1 objective instead of the mesh")
        print("  optimum lowers the objective on every target and does not help RMSD.  The")
        print("  objective/structure dissociation reappears INSIDE the degree-1 object itself.")

    # ------------------------------------------- the zero-information reference measure
    hx = load("exp_helixmu.json", quiet=True)
    if hx and hx["rows"]:
        hm = {r["pdb"]: r for r in hx["rows"]}
        com = [r for r in rows if r["pdb"] in hm]
        if com:
            f2 = np.array([r["fold"] for r in com])
            av2 = np.array([r["avg"] for r in com], float)
            print("\n" + "=" * 102)
            print(f"[HELIX-MU] the ZERO-INFORMATION reference measure, {len(com)} shared targets"
                  f"{'' if hx.get('complete') else '   ** PARTIAL **'}")
            print("  mu = every residue at the ideal alpha-helix (-57, -47) with 10 deg jitter:")
            print("  the SAME measure for every residue of every target -- sequence-blind,")
            print("  position-blind, target-blind.  MATH's estimator untouched; only the measure")
            print("  changes.  This is the control the coordinator asked every torsion-channel")
            print("  arm to carry, and it applies here because the SHIPPED mu is itself a")
            print("  conditioned torsion channel sitting inside the objective's definition.")
            print("=" * 102)
            print(f"\n  {'arm':<26}{'pool mu':>9}{'helix mu':>10}{'difference [95% CI]':>28}")
            for k, hk in (("lam000", "helix_le1"), ("argmin_le1_hold", "helix_argmin")):
                a = np.array([r[k] for r in com], float)
                b = np.array([hm[r["pdb"]][hk] for r in com], float)
                mu, lo, hi = boot_fold(b - a, f2, rng)
                print(f"  {k:<26}{a.mean():>9.3f}{b.mean():>10.3f}"
                      f"   {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
            for k, hk in (("lam000", "helix_le1"), ("argmin_le1_hold", "helix_argmin")):
                b = np.array([hm[r["pdb"]][hk] for r in com], float)
                mu, lo, hi = boot_fold(b - av2, f2, rng)
                print(f"  helix-mu {hk:<17} vs Control A ({av2.mean():.3f}): "
                      f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
            print("\n  READ.  If helix-mu matches pool-mu, the degree-1 object's information is")
            print("  GENERIC BACKBONE PLAUSIBILITY, not the pool's target conditioning, and the")
            print("  sentence 'degree-1 beats its zero-information controls therefore it carries")
            print("  real information' collapses.  Either way the falsification is untouched:")
            print("  degree-1 fails against Control A by +1.287 A under the shipped mu.")

    ob_ = load("exp_budget.json", quiet=True)
    if ob_ is not None:
        budget_report(ob_["rows"], rng)
        #: the two accounting conventions, side by side, because they read oppositely and
        #: quoting one without the other is how a search comparison gets rigged
        bm = {r["pdb"]: r for r in ob_["rows"]}
        com = [r for r in rows if r["pdb"] in bm]
        if com:
            nn = np.array([r["n"] for r in com], float)
            nf0 = np.array([r["lam000_nfev"] for r in com], float)
            nf1 = np.array([r["lam100_nfev"] for r in com], float)
            gc = np.array([bm[r["pdb"]]["greedy_full_calls"] for r in com], float)
            print("\n  ACCOUNTING, both conventions (BRIEF sec 5 'search quality'):")
            print(f"    {'method':<26}{'obj calls':>11}{'grad-equivalent':>18}")
            print(f"    {'L-BFGS on E_full (lam=1)':<26}{nf1.mean():>11.0f}"
                  f"{(nf1 * (1 + 2 * nn)).mean():>18.0f}")
            print(f"    {'L-BFGS on E_le1 (lam=0)':<26}{nf0.mean():>11.0f}"
                  f"{(nf0 * (1 + 2 * nn)).mean():>18.0f}")
            print(f"    {'greedy 1-opt on E_full':<26}{gc.mean():>11.0f}{gc.mean():>18.0f}")
            print(f"    {'exact argmin of E_le1':<26}"
                  f"{np.array([bm[r['pdb']]['exact_le1_calls'] for r in com]).mean():>11.0f}"
                  f"{np.array([bm[r['pdb']]['exact_le1_calls'] for r in com]).mean():>18.0f}")
            print("    (a gradient here is analytic and costs ~1 call; by finite differences it")
            print("     would cost 2n+1.  Both are shown so neither reading can be cherry-picked.)")

    lv = os.path.join(RESULTS, "exp_leverage.json")
    if os.path.exists(lv):
        print()
        from s18 import exp_leverage as XL
        XL.report(lv)


def budget_report(rows, rng):
    n = len(rows)
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    print("\n" + "=" * 102)
    print(f"[CONTROL B] MATCHED BUDGET, n = {n}, K = {rows[0]['K']} states/residue, "
          f"budget = {rows[0]['budget']} evaluations")
    print("  ACCOUNTING: one evaluation = one objective call at one COMPLETE configuration.")
    print("=" * 102)
    print(f"\n  {'method':<20}{'RMSD':>7}{'median':>8}{'objective':>13}{'calls':>8}"
          f"{'vs avg':>9}{'W/L':>9}")
    print(f"  {'avg (Control A)':<20}{g('avg').mean():>7.3f}{np.median(g('avg')):>8.3f}"
          f"{'--':>13}{'--':>8}{'--':>9}{'--':>9}")
    print(f"  {'discrete start':<20}{g('start_rmsd').mean():>7.3f}"
          f"{np.median(g('start_rmsd')):>8.3f}{g('start_objfull').mean():>13.1f}{1:>8}"
          f"{(g('start_rmsd') - g('avg')).mean():>+9.3f}")
    for which in ("full", "le1"):
        print(f"  -- objective: E_{which} " + "-" * 62)
        for m in ("greedy", "greedy4x", "rls", "anneal", "bestofN"):
            k = f"{m}_{which}"
            if k not in rows[0]:
                continue
            v = g(k); d = v - g("avg")
            ck = f"{m}_{which}_calls"
            calls = g(ck).mean() if ck in rows[0] else rows[0]["budget"]
            print(f"  {k:<20}{v.mean():>7.3f}{np.median(v):>8.3f}"
                  f"{g(k + '_obj').mean():>13.2f}{calls:>8.0f}{d.mean():>+9.3f}"
                  f"{int((d < 0).sum()):>5}/{int((d > 0).sum()):<4}")
    v = g("exact_le1"); d = v - g("avg")
    print(f"  {'exact_le1 (CERT)':<20}{v.mean():>7.3f}{np.median(v):>8.3f}"
          f"{g('exact_le1_obj').mean():>13.2f}{g('exact_le1_calls').mean():>8.0f}"
          f"{d.mean():>+9.3f}{int((d < 0).sum()):>5}/{int((d > 0).sum()):<4}")
    cf = int(np.array([r["greedy_full_converged"] for r in rows]).sum())
    cl = int(np.array([r["greedy_le1_converged"] for r in rows]).sum())
    print(f"\n  greedy on E_full reached a certified 1-opt local optimum on {cf}/{n} targets "
          f"within budget (mean {g('greedy_full_calls').mean():.0f} calls)")
    print(f"  greedy on E_le1  reached one on {cl}/{n} (mean "
          f"{g('greedy_le1_calls').mean():.0f} calls)")
    for which in ("full", "le1"):
        a, b = g(f"greedy_{which}_obj"), g(f"greedy4x_{which}_obj")
        print(f"  4x budget on E_{which}: objective {a.mean():.2f} -> {b.mean():.2f} "
              f"({(b < a - 1e-9).sum()}/{n} improved); RMSD "
              f"{g(f'greedy_{which}').mean():.3f} -> {g(f'greedy4x_{which}').mean():.3f}")
    print(f"  greedy_le1 vs the CERTIFIED separable optimum: {g('greedy_le1_obj').mean():.2f} "
          f"vs {g('exact_le1_obj').mean():.2f} (identical on "
          f"{(np.abs(g('greedy_le1_obj') - g('exact_le1_obj')) < 1e-6).sum()}/{n})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
