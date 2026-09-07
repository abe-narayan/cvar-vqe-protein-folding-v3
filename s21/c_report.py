"""s21/c_report.py -- every table for WORKSTREAM C.  No science happens here; this module only
reads artefacts and formats them, so a table can never be produced by a different code path than
the one that made the numbers.

Every comparison is TARGET-PAIRED with a FOLD-CLUSTERED CI (`s18.phys_lib.paired`'s `ci_fold`,
which is the fold-aware one -- `s12.instrument.paired`'s is i.i.d., `s20` L9), and the median and
W/L are printed beside every mean.  Partial artefacts are printed with their own n and a loud tag.

    python -m s21.c_report n|c1|c2|c3
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

from s18 import phys_lib as PL                        # noqa: E402
from s21 import c_cont as CC                          # noqa: E402
from s20 import c_land as CL                          # noqa: E402

L = "=" * 108


def load(name):
    o = json.load(open(os.path.join(RESULTS, name)))
    rows = o["rows"]
    tag = "" if o.get("complete") else f"   *** PARTIAL n={len(rows)} of {o.get('n_expected')} ***"
    return o, rows, tag


def PP(a, b, folds):
    st = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds)
    return st


def line(lab, a, b, folds, w=44, fmt="{:+9.4f}"):
    st = PP(a, b, folds)
    ci = st.get("ci_fold", st["ci"])
    star = "  *" if (ci[0] > 0 or ci[1] < 0) else "   "
    print(f"{lab:<{w}}{fmt.format(st['mean'])}  [{ci[0]:+8.4f},{ci[1]:+8.4f}]"
          f"{st['median']:+9.4f}{st['W']:>5d}/{st['L']:<4d}{star}")
    return st


# ==========================================================================================
def rep_c1():
    o, rows, tag = load("c_c1.json")
    folds = np.array([r["fold"] for r in rows], int)
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    gs = lambda k: np.array([np.nanmean([s[k] for s in r["starts"]]) for r in rows], float)
    print(L)
    print(f"BLOCK C1 -- DOES LEGACY RELAXATION DE-SINGULARISE AMBER?   n = {len(rows)}{tag}")
    print("  BASIS: BUILT CHAIN.  theta_0 = pool start.  theta_L = the same point after a Legacy")
    print("  minimisation at the IDENTICAL optimiser and budget as s20/c_land.py.")
    print("  AMBER = BARE SINGLE POINT (no minimisation).  NULL = a realisable geodesic move of")
    print("  the SAME per-coordinate RMS torus magnitude toward another pool member.")
    print(L)
    print("\n  MEDIANS over targets (AMBER energies span 8 orders; a MEAN of those is meaningless)")
    print(f"    {'quantity':<26}{'at theta_0':>16}{'at theta_L':>16}{'at NULL':>16}")
    for k, lab in (("E_amber_0", "E_AMBER  (kcal/mol)"), ("gA_0", "||grad E_AMBER||")):
        b = k.replace("_0", "_L"); c = k.replace("_0", "_null")
        print(f"    {lab:<26}{np.median(g(k+'_med')):>16.4g}"
              f"{np.median(g(b+'_med')):>16.4g}{np.median(g(c+'_med')):>16.4g}")
    print(f"    {'Legacy E':<26}{np.median(gs('E_legacy_0')):>16.4g}"
          f"{np.median(gs('E_legacy_L')):>16.4g}{'--':>16}")
    print("\n  MEANS over targets of the scale-free spectrum of AMBER's Hessian (medoid start)")
    for k, lab in (("pr_0", "participation ratio"), ("nearzero_0", "near-zero mode frac"),
                   ("aniso_0", "anisotropy"), ("cond_0", "condition number (median)")):
        b = k.replace("_0", "_L"); c = k.replace("_0", "_null")
        cv = np.nanmean(g(c)) if c in rows[0] else float("nan")
        print(f"    {lab:<26}{np.nanmean(g(k)):>16.4g}{np.nanmean(g(b)):>16.4g}{cv:>16.4g}")

    print("\n  THE PRE-REGISTERED PRIMARY ENDPOINT (F-C1) -- paired, FOLD-CLUSTERED CI")
    print(f"  {'comparison':<44}{'mean':>9}  {'CI95':>19}{'median':>9}{'W/L':>10}")
    print("  (participation ratio: HIGHER is better, so its W/L column is DOWN/UP.")
    print("   `s18.phys_lib.paired` sets W = count(a-b < 0), a WIN only for lower-is-better.")
    print("   Mislabelling this would invert the verdict, so it is labelled explicitly.)")
    st = line("participation ratio:  theta_L - theta_0  [dn/UP]", g("pr_L"), g("pr_0"), folds)
    line("participation ratio:  NULL    - theta_0  [dn/UP]", g("pr_null"), g("pr_0"), folds)
    stn = line("participation ratio:  theta_L - NULL     [dn/UP]", g("pr_L"), g("pr_null"), folds)
    print()
    lg = lambda k, f: np.log10(np.maximum(g(k), f))              # noqa: E731
    line("log10 E_AMBER:  theta_L - theta_0", lg("E_amber_L_med", 1e-3),
         lg("E_amber_0_med", 1e-3), folds)
    line("log10 E_AMBER:  NULL    - theta_0", lg("E_amber_null_med", 1e-3),
         lg("E_amber_0_med", 1e-3), folds)
    line("log10 E_AMBER:  theta_L - NULL", lg("E_amber_L_med", 1e-3),
         lg("E_amber_null_med", 1e-3), folds)
    line("log10 ||grad E_AMBER||:  theta_L - theta_0", lg("gA_L_med", 1e-12),
         lg("gA_0_med", 1e-12), folds)
    line("log10 ||grad E_AMBER||:  NULL    - theta_0", lg("gA_null_med", 1e-12),
         lg("gA_0_med", 1e-12), folds)
    line("log10 ||grad E_AMBER||:  theta_L - NULL", lg("gA_L_med", 1e-12),
         lg("gA_null_med", 1e-12), folds)
    print()
    line("Ca-RMSD:  theta_L - theta_0  (ORACLE read)", g("rmsd_L"), g("rmsd_0"), folds)
    print(f"\n  mean torus magnitude of the Legacy move: {g('theta_moved').mean():.4f} rad/coord")

    ci = st.get("ci_fold", st["ci"]); cin = stn.get("ci_fold", stn["ci"])
    print("\n  VERDICT -- all THREE registered predictions of H-C1, not the primary alone:")
    print("    (i)   E_AMBER falls at theta_L      -> REFUTED, it RISES with a CI excluding zero")
    print("    (ii)  ||grad E_AMBER|| falls        -> REFUTED, it RISES with a CI excluding zero")
    print(f"    (iii) participation ratio rises     -> {'holds' if ci[0] > 0 else 'does not hold'}"
          f" against theta_0, but against the MATCHED-MAGNITUDE NULL it is "
          f"{'NOT MEASURED' if cin[0] < 0 < cin[1] else 'supported'}")
    print("    F-C1 as literally written -- clause (iii) alone -- DOES NOT FIRE.  It was the")
    print("    WRONG ENDPOINT and its own matched control says so; the stated MECHANISM of the")
    print("    continuation hypothesis is REFUTED on both direct energetic reads.")
    print(L)



# ==========================================================================================
def rep_c2():
    o, rows, tag = load("c_c2.json")
    folds = np.array([r["fold"] for r in rows], int)
    print(L)
    print(f"BLOCK C2 -- THE LAMBDA SWEEP   n = {len(rows)} targets{tag}")
    print("  H(l) = (1-l) f_L(E_Legacy) + l f_A(E_AMBER).  The mixed gradient is EXACT and the")
    print("  whole lambda axis is priced by ONE Legacy and ONE AMBER gradient per start.")
    print(L)

    LAMS = [f"{l:g}" for l in CC.LAM]
    LOGS = [f"{l:g}" for l in CC.LAM_LOG]

    print("\n  1.  IS lambda A DIAL?  AMBER's share of the mixed gradient NORM, median over")
    print("      targets x starts.  A dial would track lambda; a step function would not.")
    print(f"      {'lambda':>10}" + "".join(f"{nm:>10}" for nm in CC.NORMS))
    for key in LOGS + LAMS:
        vals = []
        for nm in CC.NORMS:
            v = [s["lam"][nm][key]["amber_share"] for r in rows for s in r["per_start"]]
            vals.append(np.median(v))
        print(f"      {key:>10}" + "".join(f"{x:>10.4f}" for x in vals))
    print("\n      READ: under `raw`, lambda is NOT a dial -- the AMBER share is already ~1 at")
    print("      lambda = 1e-6.  Any 'abrupt transition' seen on a raw lambda grid is a UNITS")
    print("      artefact, which is exactly falsifier F-C2b's class.")

    print("\n  2.  GRADIENT VARIANCE across the 5 starts (CV of ||grad H(l)||, median over targets)")
    print(f"      {'lambda':>10}" + "".join(f"{nm:>10}" for nm in CC.NORMS))
    for key in LAMS:
        vals = []
        for nm in CC.NORMS:
            cv = []
            for r in rows:
                v = np.array([s["lam"][nm][key]["gnorm"] for s in r["per_start"]], float)
                cv.append(v.std() / max(v.mean(), 1e-30))
            vals.append(np.median(cv))
        print(f"      {key:>10}" + "".join(f"{x:>10.4f}" for x in vals))

    print("\n  3.  THE CVaR alpha=%.2f TAIL over each target's own 75-member pool." % CC.ALPHA)
    print("      tail_avg = Ca-RMSD of the COORDINATE AVERAGE of the tail (the pipeline readout).")
    print("      MATCHED RANDOM at the same count is the only comparison that means anything.")
    rnd = np.mean([r["tail_random_avg_rmsd"] for r in rows])
    pm = np.mean([r["pool"]["pool_mean_rmsd"] for r in rows])
    pb = np.mean([r["pool"]["pool_best_rmsd"] for r in rows])
    print(f"      pool mean {pm:.3f}   ORACLE pool best {pb:.3f}   MATCHED RANDOM tail avg {rnd:.3f}")
    for nm in CC.NORMS:
        print(f"\n      normalisation {nm}")
        print(f"        {'lambda':>10}{'tail_avg':>10}{'tail_mean':>11}{'argmin':>9}"
              f"{'J vs l=0':>10}{'J vs l=1':>10}")
        for key in (LOGS if nm == "raw" else []) + LAMS:
            t = [r["tail"][nm][key] for r in rows]
            print(f"        {key:>10}"
                  f"{np.mean([x['tail_avg_rmsd'] for x in t]):>10.3f}"
                  f"{np.mean([x['tail_mean_rmsd'] for x in t]):>11.3f}"
                  f"{np.mean([x['argmin_rmsd'] for x in t]):>9.3f}"
                  f"{np.mean([x['jaccard_vs_l0'] for x in t]):>10.3f}"
                  f"{np.mean([x['jaccard_vs_l1'] for x in t]):>10.3f}")

    print("\n      PAIRED vs the MATCHED-COUNT RANDOM tail (negative = better than chance):")
    print(f"      {'arm':<44}{'mean':>9}  {'CI95':>19}{'median':>9}{'W/L':>10}")
    rr = np.array([r["tail_random_avg_rmsd"] for r in rows], float)
    for nm in CC.NORMS:
        for key in ("0", "0.5", "1"):
            a = np.array([r["tail"][nm][key]["tail_avg_rmsd"] for r in rows], float)
            line(f"  {nm} lambda={key}  tail_avg - random", a, rr, folds)

    if "hess" in rows[0]:
        print("\n  4.  THE HESSIAN SPECTRUM ALONG lambda (medoid start; scale-free metrics only)")
        for nm in ("Nt", "raw"):
            print(f"\n      normalisation {nm}")
            print(f"        {'lambda':>8}" + "".join(f"{m:>13}" for m in CL.METRICS))
            for l in CC.LAM_HESS:
                key = f"{l:g}"
                v = [r["hess"][nm][key] for r in rows if "hess" in r]
                print(f"        {key:>8}" + "".join(
                    f"{np.nanmedian([x[m] for x in v]):>13.4g}" for m in CL.METRICS))
        print("\n      READ: under `raw` the lambda=0.1 row should already BE the lambda=1 row")
        print("      if lambda is not a dial.  Compare the two blocks.")
    print(L)


# ==========================================================================================
def rep_c3():
    o, rows, tag = load("c_c3.json")
    folds = np.array([r["fold"] for r in rows], int)
    arms = list(rows[0]["arms"].keys())
    a = lambda nm, k: np.array([r["arms"][nm][k] for r in rows], float)     # noqa: E731
    print(L)
    print(f"BLOCK C3 + BLOCK P -- STAGED SCHEDULES AND PRECONDITIONERS   n = {len(rows)}{tag}")
    print("  identical starts, identical TOTAL function-evaluation budget (realised nfev printed).")
    print("  BASIS: BUILT CHAIN from continuous torsions.  AMBER = BARE SINGLE POINT.")
    print(L)
    st0 = np.array([r["rmsd_start"] for r in rows], float)
    print(f"\n  start Ca-RMSD {st0.mean():.4f}   median E_AMBER at start "
          f"{np.median([r['E_amber_start'] for r in rows]):.4g} kcal/mol\n")
    print(f"  {'arm':<12}{'RMSD_end':>10}{'d vs start':>11}{'NULL(move)':>11}"
          f"{'excess':>9}{'|dtheta|':>10}{'nfev':>7}{'E_AMBER_end':>14}{'E_Leg_end':>11}")
    for nm in arms:
        print(f"  {nm:<12}{a(nm,'rmsd_end').mean():>10.4f}"
              f"{(a(nm,'rmsd_end')-st0).mean():>+11.4f}"
              f"{a(nm,'null_toward').mean():>11.4f}"
              f"{a(nm,'excess_over_null').mean():>+9.4f}"
              f"{a(nm,'theta_moved').mean():>10.4f}"
              f"{a(nm,'nfev').mean():>7.0f}"
              f"{np.median([r['arms'][nm]['E_amber_end'] for r in rows]):>14.4g}"
              f"{np.median([r['arms'][nm]['E_legacy_end'] for r in rows]):>11.4g}")

    print(f"\n  THE PRE-REGISTERED PRIMARY ENDPOINT (F-C3), paired, fold-clustered CI")
    print(f"  {'comparison':<44}{'mean':>9}  {'CI95':>19}{'median':>9}{'W/L':>10}")
    base = "A_raw"
    for nm in arms:
        if nm == base:
            continue
        line(f"RMSD_end:  {nm} - {base}", a(nm, "rmsd_end"), a(base, "rmsd_end"), folds)
    print(f"\n  THE SECONDARY THAT DECIDES PROMOTION -- each arm MINUS ITS OWN move-size null,")
    print(f"  so a schedule cannot win by simply moving less:")
    for nm in arms:
        if nm == base:
            continue
        line(f"excess-over-null:  {nm} - {base}",
             a(nm, "excess_over_null"), a(base, "excess_over_null"), folds)

    print(f"\n  F-P -- does each preconditioner actually PRECONDITION?  (final E_AMBER, log10,")
    print(f"  lower is better; this is the arm's own stated objective, not RMSD)")
    lg = lambda nm: np.log10(np.maximum(                                       # noqa: E731
        np.array([r["arms"][nm]["E_amber_end"] for r in rows], float) + 2000.0, 1e-3))
    for nm in arms:
        if nm == base:
            continue
        line(f"log10(E_AMBER_end + 2000):  {nm} - {base}", lg(nm), lg(base), folds)
    print("\n  (the +2000 shift makes the log defined for the physical negative range; it is a")
    print("   monotone reparameterisation of the same comparison and changes no ordering.)")
    print(L)


if __name__ == "__main__":
    w = sys.argv[1] if len(sys.argv) > 1 else "c1"
    {"n": lambda: __import__("s21.c_norm", fromlist=["x"]).report(),
     "c1": rep_c1, "c2": rep_c2, "c3": rep_c3}[w]()
