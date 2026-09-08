"""s18/phys_report.py -- the PHYSICS workstream's tables.

Every table obeys the sprint's standing law: TARGET is the unit, a mean is never printed without
its median and win/loss, every interval is the FOLD-CLUSTERED one, every arm carries both
mandatory controls, and every AMBER arm is compared to **its own gated input** with the exclusion
count and the excluded PDB IDs printed beside it.

Two substitutions the brief forbids, and this module is built so neither can happen silently:
  * an AMBER validity statistic never appears without the RMSD of the SAME structures;
  * a Legacy correlation never appears without the selection delta against a matched-random
    operation of the same count.
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

from s18 import phys_lib as PL              # noqa: E402

#: The coordinator's n = 126 reference scale (`s18/results/objceil.json`, COORD_FINDING.md).
#: A lambda arm has no interpretation without it.
REF = {"coordinate average (start)": 3.048, "alpha=0 (deployed distogram)": 3.610,
       "shuffled (direction destroyed)": 2.609, "isotropic (matched noise)": 2.573,
       "alpha=0.5": 2.448, "alpha=1 (perfect distances)": 1.152}


def _load(name):
    p = os.path.join(PL.RESULTS, name)
    with open(p) as fh:
        o = json.load(fh)
    rows = o.get("rows", [])
    if not o.get("complete", False):
        print(f"  !! {name}: PARTIAL, {len(rows)} rows -- read as a smoke, never as a result")
    return o, rows


def _g(rows, path, default=np.nan):
    """Pull a possibly-nested key from every row; missing -> nan."""
    out = []
    for r in rows:
        v = r
        try:
            for k in path:
                v = v[k]
            out.append(float(v))
        except (KeyError, TypeError, IndexError):
            out.append(default)
    return np.asarray(out, float)


def _line(nm, v, base, folds, names, extra="", seed=3):
    p = PL.paired(v, base, folds=folds, names=names, seed=seed)
    ci = p.get("ci_fold", p["ci"])
    return (f"  {nm:<26}{np.nanmean(v):>8.3f}{np.nanmedian(v):>9.3f}   "
            f"{p['mean']:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]{p['W']:>6}/{p['L']}{extra}")


# ==========================================================================
# PHASE 6 -- the lambda ladder
# ==========================================================================
def report_lambda(name="lam.json"):
    o, rows = _load(name)
    cfg = o.get("config", {})
    rows = [r for r in rows if "refine" in r]
    folds = np.array([r["fold"] for r in rows], int)
    names = [r["pdb"] for r in rows]
    print("=" * 100)
    print("PHASE 6 -- `leg_contact` AS A TERM IN THE OBJECTIVE.   n = %d" % len(rows))
    print("=" * 100)
    print(f"  degree-1 object      : {cfg.get('anova_source')}  cfg {cfg.get('anova_hash')}")
    print(f"  normalisation        : {cfg.get('normalisation')}")
    print(f"  ladder (PRE-REGISTERED, not extended): {cfg.get('LAM')}   PRIMARY lam_c = "
          f"{cfg.get('lam_primary')}")

    # ---- the normalisation, printed before any lambda number
    s = _g(rows, ("scale", "s"))
    sd_d = _g(rows, ("scale", "sd_dist"))
    sd_c = _g(rows, ("scale", "sd_contact"))
    s75 = _g(rows, ("scale_top75", "s"))
    ndeg = int(sum(1 for r in rows if r.get("scale", {}).get("degenerate")))
    print(f"\n  THE SCALE (native-free, K=500 ranker-neutral pool).  s = sd(E_le1)/sd(E_contact)")
    print(f"    sd(E_le1)     median {np.nanmedian(sd_d):.4g}   sd(E_contact) median {np.nanmedian(sd_c):.4g}")
    print(f"    s             median {np.nanmedian(s):.4g}   IQR [{np.nanpercentile(s,25):.4g}, "
          f"{np.nanpercentile(s,75):.4g}]   degenerate {ndeg}")
    print(f"    SENSITIVITY: the same scale on the shipped top-75 is median {np.nanmedian(s75):.4g}"
          f"  (ratio to pool {np.nanmedian(s75/np.where(s>0,s,np.nan)):.3f})")

    # ---- CROSS-INSTRUMENT CHECK against the coordinator's alpha = 0 arm
    try:
        oc = json.load(open(os.path.join(PL.RESULTS, "objceil.json")))
        m = {r["pdb"]: r for r in oc.get("rows", oc.get("per_target", []))}
        a0 = np.array([m[r["pdb"]]["a0.0"] if r["pdb"] in m else np.nan for r in rows], float)
        mine = _g(rows, ("refine", "full", "rmsd"))
        st = np.array([m[r["pdb"]]["proj"] if r["pdb"] in m else np.nan for r in rows], float)
        ok = np.isfinite(a0) & np.isfinite(mine)
        print(f"\n  CROSS-INSTRUMENT CHECK (gate G2b).  My `full` arm and the coordinator's "
              f"alpha = 0 arm\n  are the SAME objective, the SAME dhat and the SAME start, "
              f"optimised by different code.")
        print(f"    n = {int(ok.sum())}   mine {np.nanmean(mine[ok]):.3f}   "
              f"objceil {np.nanmean(a0[ok]):.3f}   mean |diff| {np.nanmean(np.abs(mine-a0)[ok]):.4f}"
              f"   max |diff| {np.nanmax(np.abs(mine-a0)[ok]):.4f}")
        sp = _g(rows, ("proj",))
        print(f"    start (projection): mine {np.nanmean(sp[ok]):.4f}  objceil "
              f"{np.nanmean(st[ok]):.4f}  max |diff| {np.nanmax(np.abs(sp-st)[ok]):.2e}")
    except Exception as ex:
        print(f"  (cross-instrument check unavailable: {ex})")

    # ---- ALIGNMENT
    print("\n" + "-" * 100)
    print("  6a  OBJECTIVE ALIGNMENT -- no optimisation.  K = 500 ranker-neutral pool.")
    print("      `sel` is `argmin_tied` (the tie set is AVERAGED, never read in cache order).")
    print("-" * 100)
    keys = ["full", "le1"] + [PL_k for PL_k in [f"lc{int(round(l*100)):+04d}" for l in cfg.get("LAM", PL.LAM_LADDER)]] \
        + ["contact", "null" + f"lc{int(round(PL.LAM_PRIMARY*100)):+04d}", "rand"]
    base_sel = _g(rows, ("align", "rand", "sel_band"))
    print(f"  {'objective':<26}{'rho glob':>10}{'rho band':>10}{'sel glob':>10}{'sel band':>10}"
          f"{'  band sel vs MATCHED RANDOM [95% fold CI]':>44}{'W/L':>9}")
    for k in keys:
        if not any(k in r.get("align", {}) for r in rows):
            continue
        rg = _g(rows, ("align", k, "rho_global"))
        rb = _g(rows, ("align", k, "rho_band"))
        sg = _g(rows, ("align", k, "sel_global"))
        sb = _g(rows, ("align", k, "sel_band"))
        p = PL.paired(sb, base_sel, folds=folds, names=names, seed=3)
        ci = p.get("ci_fold", p["ci"])
        tag = "  [CONTROL]" if k in ("rand",) or k.startswith("null") else ""
        print(f"  {k:<26}{np.nanmean(rg):>10.3f}{np.nanmean(rb):>10.3f}{np.nanmean(sg):>10.3f}"
              f"{np.nanmean(sb):>10.3f}      {p['mean']:+.4f} [{ci[0]:+.4f},{ci[1]:+.4f}]"
              f"{p['W']:>6}/{p['L']}{tag}")

    # ---- REFINEMENT
    print("\n" + "-" * 100)
    print("  6b  LOCAL REFINEMENT AND FINAL RMSD -- identical start (the projected coordinate")
    print("      average).  Baseline column is lam_c = 0, i.e. DEGREE-1 DISTANCE ONLY.")
    print("-" * 100)
    for k, v in REF.items():
        print(f"      reference scale (coordinator, n=126):  {k:<38} {v:.3f}")
    lam_keys = [f"lc{int(round(l*100)):+04d}" for l in cfg.get("LAM", PL.LAM_LADDER)]
    zero = f"lc{0:+04d}"
    base = _g(rows, ("refine", zero, "rmsd"))
    start = _g(rows, ("refine", "start", "rmsd"))
    avg = _g(rows, ("avg",))
    #: THE REALISED SHARE.  `lam_c` is matched on the two terms' spread over CANDIDATES, but the
    #: optimiser drives them over a completely different range.  `share` is the fraction of the
    #: combined objective's TOTAL DROP from the start that the contact term supplies.  A nominal
    #: "one pool-sd against one pool-sd" can still be 90% contact in the thing actually minimised,
    #: and if it is, the ladder never tested a SMALL contact correction.  This column is the only
    #: thing that makes lam_c interpretable and it is printed beside every arm.
    s_arr = _g(rows, ("scale", "s"))
    c0 = _g(rows, ("refine", "start", "obj_contact"))
    b0 = _g(rows, ("refine", "start", "obj_le1"))
    print(f"\n  {'arm':<26}{'RMSD':>8}{'median':>9}   {'vs lam_c = 0 [95% fold CI]':<28}{'W/L':>8}"
          f"{'objLE1':>11}{'objFULL':>11}{'objCON':>11}{'disp':>8}{'share':>8}")
    print(f"  {'coordinate average':<26}{np.nanmean(avg):>8.3f}{np.nanmedian(avg):>9.3f}"
          f"   {'(the structure the pipeline builds)':<28}")
    print(_line("start (projection)", start, base, folds, names))
    for k in lam_keys + ["full", "null" + f"lc{int(round(PL.LAM_PRIMARY*100)):+04d}", "rand_move"]:
        if not any(k in r.get("refine", {}) for r in rows):
            continue
        v = _g(rows, ("refine", k, "rmsd"))
        o1 = np.nanmean(_g(rows, ("refine", k, "obj_le1")))
        o2 = np.nanmean(_g(rows, ("refine", k, "obj_full")))
        o3 = np.nanmean(_g(rows, ("refine", k, "obj_contact")))
        dp = np.nanmean(_g(rows, ("refine", k, "disp")))
        tag = ("  [CONTROL]" if k.startswith("null") or k == "rand_move"
               else ("  <-- PRIMARY" if k == f"lc{int(round(PL.LAM_PRIMARY*100)):+04d}" else ""))
        sh = ""
        if k.startswith("lc") and k != zero:
            lamv = float(k[2:]) / 100.0
            dcon = np.abs(lamv * s_arr * (_g(rows, ("refine", k, "obj_contact")) - c0))
            dbase = np.abs(_g(rows, ("refine", k, "obj_le1")) - b0)
            sh = f"{np.nanmedian(dcon / np.maximum(dcon + dbase, 1e-12)):>8.2f}"
        print(_line(k, v, base, folds, names,
                    extra=f"{o1:>11.1f}{o2:>11.1f}{o3:>11.2f}{dp:>8.2f}{sh}{tag}"))

    # ---- THE CONTRASTS THE VERDICT TURNS ON, each stated as its own paired test
    lp = f"lc{int(round(PL.LAM_PRIMARY * 100)):+04d}"
    print("\n  KEY CONTRASTS (negative = the first arm is better)")
    for a, b, q in ((lp, zero, "PRIMARY vs degree-1 only -- THE PRE-REGISTERED TEST"),
                    (lp, "null" + lp, "PRIMARY vs its ZERO-INFORMATION control (shuffled MJ):"
                                      " does the SEQUENCE information in MJ do anything?"),
                    (lp, "rand_move", "PRIMARY vs a MATCHED-RANDOM move of its own magnitude"),
                    ("full", zero, "full distance vs degree-1 only"),
                    ("full", "start", "full distance vs the START (L30's question, reproduced)"),
                    (zero, "start", "degree-1 only vs the START")):
        va, vb = _g(rows, ("refine", a, "rmsd")), _g(rows, ("refine", b, "rmsd"))
        if np.isnan(va).all() or np.isnan(vb).all():
            continue
        p = PL.paired(va, vb, folds=folds, names=names, seed=3)
        ci = p.get("ci_fold", p["ci"])
        print(f"    {a:<12} vs {b:<12} {p['mean_a']:>7.3f} vs {p['mean_b']:<7.3f} "
              f"{p['mean']:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]  med {p['median']:+.3f} "
              f"{p['W']:>4}W/{p['L']}L   {q}")

    # ---- GENERATION
    if any("gen" in r for r in rows):
        print("\n" + "-" * 100)
        print("  6c  CANDIDATE GENERATION, COVERAGE AND ENSEMBLE DIVERSITY.")
        print("      m = %d starts, the first top-75 windows in RETRIEVAL order (native-free," %
              int(np.nanmean(_g(rows, ("gen", "starts", "m")))))
        print("      identical across arms).  `ens` = RMSD of the ensemble's coordinate average.")
        print("-" * 100)
        gk = ["starts", "le1", f"lc{int(round(PL.LAM_PRIMARY*100)):+04d}", "full",
              "null" + f"lc{int(round(PL.LAM_PRIMARY*100)):+04d}"]
        gbase = _g(rows, ("gen", "le1", "ens_rmsd"))
        print(f"  {'arm':<26}{'set mean':>10}{'set best':>10}{'n<2A':>7}{'diversity':>11}"
              f"{'ens RMSD':>10}{'median':>9}   {'ens vs le1 [fold CI]':<26}{'W/L':>8}")
        for k in gk:
            if not any(k in r.get("gen", {}) for r in rows):
                continue
            sm = np.nanmean(_g(rows, ("gen", k, "set_mean")))
            sb = np.nanmean(_g(rows, ("gen", k, "set_best")))
            nn = np.nanmean(_g(rows, ("gen", k, "n_near")))
            dv = np.nanmean(_g(rows, ("gen", k, "diversity")))
            en = _g(rows, ("gen", k, "ens_rmsd"))
            p = PL.paired(en, gbase, folds=folds, names=names, seed=3)
            ci = p.get("ci_fold", p["ci"])
            tag = "  [CONTROL]" if k.startswith("null") else ""
            print(f"  {k:<26}{sm:>10.3f}{sb:>10.3f}{nn:>7.2f}{dv:>11.3f}{np.nanmean(en):>10.3f}"
                  f"{np.nanmedian(en):>9.3f}   {p['mean']:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]"
                  f"{p['W']:>5}/{p['L']}{tag}")


def report_bases(name="lam_bases.json"):
    """The DECLARED EXTENSION: the same contact term on the DEPLOYED objective and on MATH's
    residue-additive object, after MATH corrected which object the brief's bridge derives and
    after the coordinator showed the deployed functional form is sound."""
    o, rows = _load(name)
    rows = [r for r in rows if "arms" in r and r["arms"]]
    folds = np.array([r["fold"] for r in rows], int)
    names = [r["pdb"] for r in rows]
    lp = f"lc{int(round(PL.LAM_PRIMARY * 100)):+04d}"
    print("\n" + "=" * 100)
    print("PHASE 6 (DECLARED EXTENSION) -- `leg_contact` ON THE DEPLOYED OBJECTIVE AND ON THE")
    print("RESIDUE-ADDITIVE ANOVA.  n = %d.  Same term, same sign, same native-free pool-sd" % len(rows))
    print("normalisation recomputed against each base's own pool sd; ladder cut to {0, +1}.")
    print("=" * 100)
    for bn, lbl in (("full", "E_full  (the DEPLOYED objective)"),
                    ("res", "E_res   (MATH's residue-additive ANOVA -- the BRIEF's object)")):
        b = _g(rows, ("arms", f"{bn}_lc+000", "rmsd"))
        s_arr = _g(rows, ("scale", bn, "s"))
        print(f"\n  base = {lbl}")
        print(f"    pool scale s median {np.nanmedian(s_arr):.4g}")
        print(f"  {'arm':<26}{'RMSD':>8}{'median':>9}   {'vs lam_c = 0 [95% fold CI]':<28}"
              f"{'W/L':>8}{'objBASE':>11}{'objCON':>11}{'disp':>8}{'share':>8}")
        c0 = _g(rows, ("arms", f"{bn}_lc+000", "obj_contact"))
        b0 = _g(rows, ("arms", f"{bn}_lc+000", "obj_full" if bn == "full" else "obj_res"))
        for k, lam in ((f"{bn}_lc+000", 0.0), (f"{bn}_{lp}", PL.LAM_PRIMARY),
                       (f"{bn}_null{lp}", PL.LAM_PRIMARY)):
            v = _g(rows, ("arms", k, "rmsd"))
            if np.isnan(v).all():
                continue
            ob_ = np.nanmean(_g(rows, ("arms", k, "obj_full" if bn == "full" else "obj_res")))
            oc = np.nanmean(_g(rows, ("arms", k, "obj_contact")))
            dp = np.nanmean(_g(rows, ("arms", k, "disp")))
            sh = ""
            if lam != 0.0:
                dcon = np.abs(lam * s_arr * (_g(rows, ("arms", k, "obj_contact")) - c0))
                dbase = np.abs(_g(rows, ("arms", k, "obj_full" if bn == "full" else "obj_res")) - b0)
                sh = f"{np.nanmedian(dcon / np.maximum(dcon + dbase, 1e-12)):>8.2f}"
            tag = "  [CONTROL]" if "null" in k else ("  <-- PRIMARY" if lam else "")
            print(_line(k, v, b, folds, names, extra=f"{ob_:>11.1f}{oc:>11.2f}{dp:>8.2f}{sh}{tag}"))
    #: and the two bases against each other, which is MATH's and EXPERIMENT's question, not mine
    a = _g(rows, ("arms", "res_lc+000", "rmsd"))
    b = _g(rows, ("arms", "full_lc+000", "rmsd"))
    p = PL.paired(a, b, folds=folds, names=names, seed=3)
    ci = p.get("ci_fold", p["ci"])
    print(f"\n  (for MATH/EXPERIMENT, not a PHYSICS claim)  E_res vs E_full, same start, same "
          f"optimiser:\n    {np.nanmean(a):.3f} vs {np.nanmean(b):.3f}   "
          f"{p['mean']:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]  {p['W']}W/{p['L']}L")


# ==========================================================================
# PHASE 8 -- the causal downstream comparison
# ==========================================================================
VKEYS = (("rama_favoured", "ramaFav", 3), ("rama_outlier", "ramaOut", 3),
         ("n_clash_2A", "clash<2", 2), ("n_clash_2p6A", "clash<2.6", 2),
         ("min_heavy", "minHeavy", 3), ("cis_frac", "cis", 3), ("omega_dev", "wdev", 1))


def report_down(name="down.json", f=0.50):
    o, rows = _load(name)
    folds = np.array([r["fold"] for r in rows], int)
    names = [r["pdb"] for r in rows]
    tag = f"@{f:.2f}"
    arms = ["none", "leg_torsion" + tag, "leg_contact" + tag, "legacy" + tag,
            "amber_sp" + tag, "rand" + tag, "helix" + tag]
    print("\n" + "=" * 100)
    print("PHASE 8 -- THE CAUSAL DOWNSTREAM COMPARISON, on IDENTICAL structures.  n = %d" % len(rows))
    print("=" * 100)
    print("  Every filter sees the SAME shipped top-75 windows and drops the worst %.0f%% by its"
          % (100 * f))
    print("  own score; the survivors go through the SAME averaging operator.  `rand` is the")
    print("  matched-random gate of the same COUNT (3 stable_rng draws); `helix` is the")
    print("  zero-information score (Ca-RMSD to a constant ideal alpha-helix).")

    base = _g(rows, ("arms", "none", "avg_rmsd"))
    print("\n  -- ENSEMBLE / OUTPUT RMSD (the answer the pipeline emits) ------------------------")
    print(f"  {'filter':<26}{'avg':>8}{'median':>9}   {'vs NO FILTER [95% fold CI]':<28}{'W/L':>8}"
          f"{'proj':>9}{'m':>6}{'setMean':>9}{'setBest':>9}{'div':>8}")
    for a in arms:
        if not any(a in r["arms"] for r in rows):
            continue
        v = _g(rows, ("arms", a, "avg_rmsd"))
        pj = np.nanmean(_g(rows, ("arms", a, "proj_rmsd")))
        m = np.nanmean(_g(rows, ("arms", a, "m")))
        sm = np.nanmean(_g(rows, ("arms", a, "set_mean")))
        sb = np.nanmean(_g(rows, ("arms", a, "set_best")))
        dv = np.nanmean(_g(rows, ("arms", a, "diversity")))
        note = "  [CONTROL]" if a.startswith(("rand", "helix")) else ""
        print(_line(a, v, base, folds, names,
                    extra=f"{pj:>9.3f}{m:>6.0f}{sm:>9.3f}{sb:>9.3f}{dv:>8.3f}{note}"))

    print("\n  -- the same arms on the PROJECTED output (the incumbent's emitted structure) ------")
    bp = _g(rows, ("arms", "none", "proj_rmsd"))
    for a in arms:
        v = _g(rows, ("arms", a, "proj_rmsd"))
        if np.isnan(v).all():
            continue
        print(_line(a, v, bp, folds, names,
                    extra="  [CONTROL]" if a.startswith(("rand", "helix")) else ""))

    # ---- NEAR-NATIVE RECALL, on ONE denominator
    has = _g(rows, ("n_near_pool",)) > 0
    print(f"\n  -- NEAR-NATIVE RECALL (sub-2.0 A members kept), n = {int(has.sum())} targets that "
          f"HAVE one; ONE denominator")
    rb = np.array([r["arms"]["rand" + tag]["n_near"] / max(r["n_near_pool"], 1)
                   for r in rows], float)
    print(f"  {'filter':<26}{'recall':>9}{'random':>9}   {'recall - random [fold CI]':<28}{'W/L':>8}")
    for a in arms:
        if a == "none" or not any(a in r["arms"] for r in rows):
            continue
        v = np.array([r["arms"][a]["n_near"] / max(r["n_near_pool"], 1) for r in rows], float)
        p = PL.paired(v[has], rb[has], folds=folds[has], names=[n for n, h in zip(names, has) if h],
                      seed=5)
        ci = p.get("ci_fold", p["ci"])
        print(f"  {a:<26}{np.nanmean(v[has]):>9.3f}{np.nanmean(rb[has]):>9.3f}   "
              f"{p['mean']:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]{p['W']:>6}/{p['L']}")

    # ---- AMBER
    print("\n  -- AMBER RESTRAINED REPAIR, k = 30, N/CA/C (the incumbent protocol) ---------------")
    print("     THE CONVERGENCE GATE was declared before use: final energy finite and <= 1000")
    print("     kcal/mol.  Every arm is compared to ITS OWN GATED INPUT, never to the ungated")
    print("     instrument-wide input.  RMSD and validity are printed together, on the SAME")
    print("     structures, because neither is admissible alone.")
    amb = [a for a in arms if any("amb30" in r["arms"].get(a, {}) for r in rows)]
    print(f"\n  {'arm':<24}{'excl':>6}{'RMSD':>8}{'input':>8}{'  vs its OWN gated input [fold CI]':<34}"
          f"{'W/L':>8}{'CaDisp':>8}")
    for a in amb:
        v = _g(rows, ("arms", a, "amb30", "rmsd"))
        iv = _g(rows, ("arms", a, "amb30", "input_rmsd"))
        en = _g(rows, ("arms", a, "amb30", "energy"))
        cv = np.array([bool(r["arms"].get(a, {}).get("amb30", {}).get("converged", False))
                       for r in rows], bool)
        dd = _g(rows, ("arms", a, "amb30", "ca_disp"))
        g = PL.gated(v, iv, en, converged=cv, folds=folds, names=names, seed=9)
        p = g["gated"]
        ci = p.get("ci_fold", p["ci"])
        print(f"  {a:<24}{g['n_excluded']:>6}{p['mean_a']:>8.3f}{p['mean_b']:>8.3f}   "
              f"{p['mean']:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]{p['W']:>8}/{p['L']}"
              f"{np.nanmean(dd[cv]):>8.3f}")
        if g["n_excluded"]:
            print(f"      excluded: {' '.join(g['excluded'])}")

    print("\n  -- VALIDITY on the SAME structures, each AMBER arm above ITS OWN GATED INPUT -----")
    hdr = "".join(f"{lbl:>11}" for _k, lbl, _d in VKEYS)
    print(f"  {'structure':<30}{hdr}")
    for a in amb:
        cv = np.array([bool(r["arms"].get(a, {}).get("amb30", {}).get("converged", False))
                       for r in rows], bool)
        en = _g(rows, ("arms", a, "amb30", "energy"))
        ok = cv & np.isfinite(en) & (en <= 1000.0)
        for lbl, sub in ((a + "  (amb30)", ("valid",)), ("  ^ its own gated input", ("input_valid",))):
            cells = ""
            for k, _l, dg in VKEYS:
                vv = _g(rows, ("arms", a, "amb30") + sub + (k,))
                cells += f"{np.nanmean(vv[ok]):>11.{dg}f}"
            print(f"  {lbl:<30}{cells}")
    print("\n  -- VALIDITY of the un-repaired coordinate average of each filter's survivors -----")
    print(f"  {'structure':<30}{hdr}{'  spacing':>10}")
    for a in arms:
        if not any(a in r["arms"] for r in rows):
            continue
        cells = "".join(f"{np.nanmean(_g(rows, ('arms', a, 'avg_valid', k))):>11.{dg}f}"
                        for k, _l, dg in VKEYS)
        sp = np.nanmean(_g(rows, ("arms", a, "avg_ca_spacing")))
        print(f"  {a + '  (avg)':<30}{cells}{sp:>10.3f}")


# ==========================================================================
# PHASE 10 -- spacing restoration
# ==========================================================================
def report_space(name="space.json"):
    o, rows = _load(name)
    folds = np.array([r["fold"] for r in rows], int)
    names = [r["pdb"] for r in rows]
    print("\n" + "=" * 100)
    print("PHASE 10 (SECONDARY) -- CHEAPER SPACING RESTORATION.  n = %d" % len(rows))
    print("=" * 100)
    print("  The averaging operator contracts Ca-Ca spacing to 2.949 A against an ideal 3.80 A")
    print("  (L28), Spearman(spacing, cis) = -0.949, and the incumbent projection costs +0.164 A")
    print("  to fix it.  The question is whether anything cheaper restores the spacing.")
    base = _g(rows, ("arms", "none", "rmsd"))
    order = ["none", "scale_uniform", "bond_norm", "segment", "minproj", "torsion_rebuild",
             "proj", "proj_fit", "helix"]
    print(f"\n  {'arm':<26}{'RMSD':>8}{'median':>9}   {'vs the AVERAGE [fold CI]':<28}{'W/L':>8}"
          f"{'spacing':>9}{'disp':>8}{'  vs its MATCHED RANDOM':>24}")
    for a in order:
        if not any(a in r["arms"] for r in rows):
            continue
        v = _g(rows, ("arms", a, "rmsd"))
        sp = np.nanmean(_g(rows, ("arms", a, "spacing")))
        dp = np.nanmean(_g(rows, ("arms", a, "disp")))
        rk = "rand@" + a
        ex = ""
        if any(rk in r["arms"] for r in rows):
            rv = _g(rows, ("arms", rk, "rmsd"))
            p2 = PL.paired(v, rv, folds=folds, names=names, seed=11)
            c2 = p2.get("ci_fold", p2["ci"])
            ex = f"   {p2['mean']:+.3f} [{c2[0]:+.3f},{c2[1]:+.3f}]"
        note = "  [ZERO-INFO CONTROL]" if a == "helix" else ""
        print(_line(a, v, base, folds, names, extra=f"{sp:>9.3f}{dp:>8.3f}{ex}{note}", seed=11))
    #: the incumbent projection is the thing to beat: repeat every arm against IT.
    pb = _g(rows, ("arms", "proj", "rmsd"))
    print(f"\n  -- against the INCUMBENT PROJECTION (the +0.164 A reference) ----------------------")
    for a in order:
        if a in ("proj",) or not any(a in r["arms"] for r in rows):
            continue
        v = _g(rows, ("arms", a, "rmsd"))
        print(_line(a, v, pb, folds, names, seed=11))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "lam"):
        try:
            report_lambda()
        except FileNotFoundError:
            print("  (lam.json not present yet)")
    if which in ("all", "bases"):
        try:
            report_bases()
        except FileNotFoundError:
            print("  (lam_bases.json not present yet)")
    if which in ("all", "down"):
        try:
            report_down(f=0.50)
            print("\n\n  ###### SHAPE ONLY: the same filters at f = 0.25 (no projection, no"
                  " AMBER).\n  ###### The PRE-REGISTERED PRIMARY is f = 0.50 above and nothing"
                  " here may be promoted over it. ######")
            report_down(f=0.25)
        except FileNotFoundError:
            print("  (down.json not present yet)")
    if which in ("all", "space"):
        try:
            report_space()
        except FileNotFoundError:
            print("  (space.json not present yet)")
