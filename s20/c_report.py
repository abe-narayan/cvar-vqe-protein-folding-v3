"""s20/c_report.py -- every Sprint-20 Agent-C table, built from the persisted artefacts only.

No table here recomputes anything: each reads a `complete` artefact and prints. TARGET is the
unit everywhere; every mean is printed with its fold-aware paired bootstrap CI, its median and
its W/L, per `s20/BRIEF.md` section 8.

    python -m s20.c_report q1
    python -m s20.c_report pareto
    python -m s20.c_report land
"""
from __future__ import annotations

import os
import sys
import json

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
from s12 import instrument as I     # noqa: E402
from s18 import phys_lib as PL      # noqa: E402


# ---------------------------------------------------------------------------------------
# A SELF-CAUGHT LABELLING DEFECT.  `s12.instrument.paired`'s `ci95` is a PLAIN i.i.d.
# target-level bootstrap; its `folds` argument only adds a per-fold mean breakdown and does
# NOT cluster the resample.  Sprint 19's numbers were made with `s18.phys_lib.paired`, which
# returns BOTH an i.i.d. `ci` and a FOLD-CLUSTERED `ci_fold`, and "fold-aware" in this
# programme means the latter (BRIEF section 8).  Every CI below is therefore produced by
# `PL.paired` and the FOLD-CLUSTERED interval is the one quoted, with the i.i.d. one printed
# beside it wherever both matter.  Caught by reading `s12/instrument.py:188` rather than
# trusting the parameter name.
def PP(a, b, folds=None):
    st = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds)
    return {"mean_diff": st["mean"], "ci95": st.get("ci_fold", st["ci"]),
            "ci_iid": st["ci"], "median_diff": st["median"],
            "n_better": st["W"], "n_worse": st["L"],
            "mean_a": st["mean_a"], "mean_b": st["mean_b"], "n": st["n"]}


def load(name, need=None):
    o = json.load(open(os.path.join(RESULTS, name)))
    if need is not None and not o.get("complete"):
        raise SystemExit(f"{name} is not complete ({o.get('n_rows')}/{o.get('n_expected')})")
    return o


def cmp_(a, b, folds, name=""):
    """paired a - b, fold-aware bootstrap, median and W/L beside the mean."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    st = PP(a, b, folds=np.asarray(folds, int))
    return (f"{name:<22s} {st['mean_a']:7.4f} vs {st['mean_b']:7.4f}  "
            f"d={st['mean_diff']:+7.4f} [{st['ci95'][0]:+7.4f},{st['ci95'][1]:+7.4f}]  "
            f"med={st['median_diff']:+7.4f}  {st['n_better']:>3d}W/{st['n_worse']:<3d}L")


def q1():
    o = load("c_q1.json", need=126)
    rows = o["rows"]
    folds = np.array([r["fold"] for r in rows], int)
    print("=" * 100)
    print("Q1  DO LEGACY AND AMBER CONTAIN COMPLEMENTARY INFORMATION?   n =", len(rows))
    print("    identical candidates: the shipped top-75 ideal-geometry rebuilds, 126 targets")
    print("=" * 100)

    print("\n1.  AGREEMENT BETWEEN THE TWO POTENTIALS  (per-target Spearman over 75 candidates)")
    print(f"{'pair':<26s}{'mean':>8s}{'median':>9s}{'sd':>8s}{'q10':>8s}{'q90':>8s}"
          f"{'>0':>6s}{'|rho|>0.8':>11s}")
    for key, lab in (("rho_leg_amb", "Legacy vs AMBER"),
                     ("rho_leg_d", "Legacy vs ORACLE d"),
                     ("rho_amb_d", "AMBER vs ORACLE d"),
                     ("rho_dis_d", "distogram vs ORACLE d"),
                     ("rho_leg_dis", "Legacy vs distogram"),
                     ("rho_amb_dis", "AMBER vs distogram")):
        v = np.array([r[key] for r in rows], float)
        print(f"{lab:<26s}{v.mean():8.3f}{np.median(v):9.3f}{v.std():8.3f}"
              f"{np.percentile(v,10):8.3f}{np.percentile(v,90):8.3f}"
              f"{(v>0).mean():6.2f}{(np.abs(v)>0.8).mean():11.2f}")
    v = np.array([r["rho_leg_amb"] for r in rows], float)
    st = PP(v, np.zeros_like(v), folds=folds)
    print(f"\n    rho(Legacy, AMBER) vs 0:  {st['mean_diff']:+.4f} "
          f"[{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}]  median {np.median(v):+.4f}")
    print(f"    F-C1(a) fires iff median rho > 0.80  ->  median = {np.median(v):.3f}  "
          f"-> {'FIRES' if np.median(v) > 0.80 else 'does NOT fire'}")

    print("\n2.  CELL SIZES  (each model 'prefers' its best 25 of 75; independence -> 8.3/16.7/16.7/33.3)")
    for c in ("both", "leg_only", "amb_only", "neither"):
        s = np.array([r["cell_sizes"][c] for r in rows], float)
        print(f"    {c:<10s} mean {s.mean():6.2f}   median {np.median(s):5.1f}   "
              f"range {s.min():.0f}-{s.max():.0f}")

    print("\n3.  WHAT THE FOUR PARTITIONS LOOK LIKE  (per-target cell means, averaged over targets)")
    axes = o["config"]["AXES"]
    print(f"{'axis':<16s}" + "".join(f"{c:>13s}" for c in
                                     ("both", "leg_only", "amb_only", "neither")))
    for a in axes:
        vals = [np.nanmean([r["cell_mean"][c][a] for r in rows]) for c in
                ("both", "leg_only", "amb_only", "neither")]
        print(f"{a:<16s}" + "".join(f"{v:13.4f}" for v in vals))

    print("\n4.  THE DISAGREEMENT CONTRAST  Legacy-only MINUS AMBER-only, "
          "against a MATCHED-RANDOM partition of identical cell sizes")
    print("    (200 random partitions per target; 'excess' = observed - null mean; z = per-target"
          " standardised, then averaged)")
    print(f"{'axis':<16s}{'observed':>11s}{'null mean':>11s}{'excess':>10s}"
          f"{'CI95 (excess)':>26s}{'med':>9s}{'W/L':>10s}{'mean z':>9s}")
    fired = []
    for a in axes:
        obs = np.array([r["contrast_leg_minus_amb"][a] for r in rows], float)
        nul = np.array([r["contrast_null"][a]["mean"] for r in rows], float)
        z = np.array([r["contrast_z"][a] for r in rows], float)
        st = PP(obs, nul, folds=folds)
        excl = (st["ci95"][0] > 0) or (st["ci95"][1] < 0)
        if excl:
            fired.append(a)
        print(f"{a:<16s}{np.nanmean(obs):11.4f}{np.nanmean(nul):11.4f}"
              f"{st['mean_diff']:10.4f}  [{st['ci95'][0]:+9.4f},{st['ci95'][1]:+9.4f}]"
              f"{st['median_diff']:9.4f}{st['n_better']:>5d}/{st['n_worse']:<4d}"
              f"{np.nanmean(z):9.2f}"
              + ("  *" if excl else ""))
    print(f"\n    axes whose Legacy-only - AMBER-only contrast has a CI EXCLUDING ZERO: "
          f"{fired if fired else 'NONE'}")
    print(f"    F-C1(b) fires iff NO axis separates -> "
          f"{'FIRES' if not fired else 'does NOT fire'}")
    print("    NOTE the axes constant by construction on an ideal-geometry rebuild "
          "(bond_strain, angle_strain, cis_frac, chirality_L_frac) carry no information here "
          "and are NOT reported as nulls.")

    print("\n4b. STRATIFICATION (BRIEF section 8: never a mean alone).  The two headline "
          "quantities, by LENGTH tercile and by FOLD.")
    ns = np.array([r["n"] for r in rows], int)
    rho = np.array([r["rho_leg_amb"] for r in rows], float)
    rgo = np.array([r["contrast_leg_minus_amb"]["rg"] for r in rows], float)
    rgn = np.array([r["contrast_null"]["rg"]["mean"] for r in rows], float)
    q1_, q2_ = np.percentile(ns, [33.3, 66.7])
    strata = [("len <= %d" % q1_, ns <= q1_),
              ("len %d-%d" % (q1_, q2_), (ns > q1_) & (ns <= q2_)),
              ("len > %d" % q2_, ns > q2_)]
    strata += [(f"fold {f}", folds == f) for f in sorted(set(folds.tolist()))]
    print(f"{'stratum':<14s}{'n':>5s}{'rho(Leg,AMB)':>14s}{'rg contrast':>14s}"
          f"{'CI95 (rg excess)':>26s}")
    for lab, m in strata:
        if m.sum() < 5:
            continue
        st = PP(rgo[m], rgn[m], folds=folds[m])
        print(f"{lab:<14s}{int(m.sum()):>5d}{rho[m].mean():+14.4f}{st['mean_diff']:+14.4f}"
              f"  [{st['ci95'][0]:+10.4f},{st['ci95'][1]:+10.4f}]")

    print("\n5.  THE MANDATED ABLATION TABLE (directive section 23)")
    print("    BASIS: POINT CLOUD -- the readout is the deployed coordinate average of the m=38")
    print("    survivors.  Comparable to Sprint 19's 3.050 A column and to nothing else.")
    print("    Validity is measured on the MEMBERS (built chains), never on the contracted mean.")
    arms = ["none", "rand", "rand2", "disto", "legacy", "amber", "leg_then_amb", "amb_then_leg"]
    keys = ["rmsd", "best_member", "mean_member", "D", "mean_align",
            "v_min_heavy", "v_n_clash_2A", "v_n_clash_2p6A",
            "v_rama_favoured", "v_rama_outlier"]
    hdr = {"rmsd": "readout", "best_member": "best mem", "mean_member": "mean mem", "D": "D",
           "mean_align": "ALIGN", "v_min_heavy": "min_heavy", "v_n_clash_2A": "clash2.0",
           "v_n_clash_2p6A": "clash2.6", "v_rama_favoured": "rama_fav",
           "v_rama_outlier": "rama_out"}
    print(f"{'arm':<15s}" + "".join(f"{hdr[k]:>11s}" for k in keys))
    A = {a: {k: np.array([r["ablation"][a][k] for r in rows], float) for k in keys} for a in arms}
    for a in arms:
        print(f"{a:<15s}" + "".join(f"{A[a][k].mean():11.4f}" for k in keys))
    print("\n    cis_frac / bond_strain / angle_strain / chirality_L_frac are EXACTLY constant "
          "across every arm (ideal-geometry rebuild); omitted rather than printed as zeros.")

    print("\n6.  EVERY ORDERED ARM AGAINST ITS MATCHED-RANDOM CONTROL  (readout Ca-RMSD, "
          "point-cloud basis)")
    for a, ctrl in (("disto", "rand"), ("legacy", "rand"), ("amber", "rand"),
                    ("leg_then_amb", "rand2"), ("amb_then_leg", "rand2")):
        print("    " + cmp_(A[a]["rmsd"], A[ctrl]["rmsd"], folds, f"{a} - {ctrl}"))
    print("\n    and against NO gate:")
    for a in arms[1:]:
        print("    " + cmp_(A[a]["rmsd"], A["none"]["rmsd"], folds, f"{a} - none"))
    print("\n    composites against their OWN FIRST STAGE:")
    print("    " + cmp_(A["leg_then_amb"]["rmsd"], A["legacy"]["rmsd"], folds,
                        "leg_then_amb - legacy"))
    print("    " + cmp_(A["amb_then_leg"]["rmsd"], A["amber"]["rmsd"], folds,
                        "amb_then_leg - amber"))
    print("=" * 100)


def _valid_cols():
    return ("min_heavy", "n_clash_2A", "n_clash_2p6A", "bond_strain", "angle_strain",
            "rama_favoured", "rama_allowed", "rama_outlier", "cis_frac", "chirality_L_frac",
            "omega_dev", "geom_rms_rel_dev")


def pareto(fn="c_pareto.json"):
    o = json.load(open(os.path.join(RESULTS, fn)))
    rows = o["rows"]
    n = len(rows)
    tag = "" if o.get("complete") else f"  *** PARTIAL n={n}/{o.get('n_expected')} ***"
    folds = np.array([r["fold"] for r in rows], int)
    arms = list(rows[0]["arms"].keys())
    print("=" * 118)
    print(f"Q3  THE AMBER REPAIR PARETO   n = {n}{tag}")
    print("    BASIS: REPAIRED EMISSION.  input = the coordinate average (point cloud); output =")
    print("    the AMBER-relaxed all-atom structure.  Every delta is against the arm's OWN input.")
    print("=" * 118)
    inp = np.array([r["input_rmsd"] for r in rows], float)
    iv = {c: np.array([r["input_valid"].get(c, np.nan) for r in rows], float)
          for c in _valid_cols()}
    print(f"\ninput (the coordinate average, POINT CLOUD): Ca-RMSD {inp.mean():.4f}")
    print("    input validity: " + "  ".join(f"{c}={np.nanmean(iv[c]):.4f}"
                                              for c in _valid_cols()))

    print(f"\n{'arm':<20s}{'Ca-RMSD':>9s}{'d vs input':>22s}{'ca_disp':>9s}"
          f"{'minheavy':>9s}{'cl2.0':>7s}{'cl2.6':>7s}{'bond':>8s}{'angle':>8s}"
          f"{'ramaF':>7s}{'ramaO':>7s}{'cis':>7s}{'conv':>6s}{'wall':>7s}")
    out = {}
    for a in arms:
        r_ = np.array([r["arms"][a]["rmsd"] for r in rows], float)
        st = PP(r_, inp, folds=folds)
        v = {c: np.array([r["arms"][a]["valid"].get(c, np.nan) for r in rows], float)
             for c in _valid_cols()}
        cd = np.array([r["arms"][a]["ca_disp"] for r in rows], float)
        cv = np.array([bool(r["arms"][a]["converged"]) for r in rows])
        wl = np.array([r["arms"][a].get("wall", 0.0) for r in rows], float)
        out[a] = {"rmsd": r_, "valid": v, "conv": cv, "disp": cd}
        print(f"{a:<20s}{r_.mean():9.4f}"
              f"  {st['mean_diff']:+6.4f}[{st['ci95'][0]:+6.4f},{st['ci95'][1]:+6.4f}]"
              f"{cd.mean():9.4f}{v['min_heavy'].mean():9.4f}"
              f"{v['n_clash_2A'].mean():7.3f}{v['n_clash_2p6A'].mean():7.2f}"
              f"{v['bond_strain'].mean():8.4f}{v['angle_strain'].mean():8.4f}"
              f"{v['rama_favoured'].mean():7.3f}{v['rama_outlier'].mean():7.3f}"
              f"{v['cis_frac'].mean():7.4f}{cv.mean():6.2f}{wl.mean():7.1f}")

    print("\nCONVERGENCE GATE (declared in PREREG_C.md section 4 before use): "
          "`core.amber.convergence_flags`")
    for a in arms:
        cv = out[a]["conv"]
        bad = [rows[i]["pdb"] for i in range(n) if not cv[i]]
        if len(bad):
            print(f"    {a:<20s} EXCLUSIONS {len(bad):>3d}/{n}: {bad[:8]}"
                  f"{' ...' if len(bad) > 8 else ''}")
    allbad = sorted({rows[i]["pdb"] for a in arms for i in range(n)
                     if not out[a]["conv"][i]})
    print(f"    union of non-converged targets over all arms: {len(allbad)}  {allbad[:12]}")
    if allbad:
        keep = np.array([r["pdb"] not in allbad for r in rows])
        print(f"\n  the same table with those {len(allbad)} targets EXCLUDED (n={keep.sum()}):")
        print(f"  {'arm':<20s}{'Ca-RMSD':>9s}{'d vs input':>22s}")
        for a in arms:
            st = PP(out[a]["rmsd"][keep], inp[keep], folds=folds[keep])
            print(f"  {a:<20s}{out[a]['rmsd'][keep].mean():9.4f}"
                  f"  {st['mean_diff']:+6.4f}[{st['ci95'][0]:+6.4f},{st['ci95'][1]:+6.4f}]")

    print("\nAGAINST THE INCUMBENT k30, on BOTH axes.  An arm is a Pareto point only if it")
    print("improves Ca-RMSD and degrades NO axis of the validity vector.")
    k30 = out["k30"]
    print(f"{'arm':<20s}{'d Ca-RMSD vs k30':>34s}{'  validity axes DEGRADED vs k30':<40s}")
    for a in arms:
        if a == "k30":
            continue
        st = PP(out[a]["rmsd"], k30["rmsd"], folds=folds)
        deg = []
        for c in _valid_cols():
            x, y = out[a]["valid"][c], k30["valid"][c]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 10 or np.allclose(x[ok], y[ok]):
                continue
            s2 = PP(x[ok], y[ok], folds=folds[ok])
            worse_is_high = c not in ("min_heavy", "rama_favoured")
            if s2["ci95"][0] > 0 or s2["ci95"][1] < 0:
                bad = (s2["mean_diff"] > 0) == worse_is_high
                if bad:
                    deg.append(f"{c}{s2['mean_diff']:+.3g}")
        print(f"{a:<20s}{st['mean_diff']:+10.4f}"
              f"[{st['ci95'][0]:+7.4f},{st['ci95'][1]:+7.4f}] med{st['median_diff']:+7.4f} "
              f"{st['n_better']:>3d}W/{st['n_worse']:<3d}L  "
              + (", ".join(deg) if deg else "NONE"))

    fnul = o.get("frame_null")
    if fnul:
        print(f"\nROTATED-FRAME NULL (exactly zero by construction; reported with its MAXIMUM): "
              f"n={fnul['n']}  MAX |d Ca-RMSD| = {fnul['max_abs']:.3e} A  "
              f"MAX |d E| = {fnul['max_abs_energy']:.3e} kcal/mol   (mean_abs "
              f"{fnul['mean_abs']:.3e} A, printed only beside the maximum)")
    print("=" * 118)


def land(fn="c_land.json"):
    o = json.load(open(os.path.join(RESULTS, fn)))
    rows = o["rows"]
    n = len(rows)
    tag = "" if o.get("complete") else f"  *** PARTIAL n={n}/{o.get('n_expected')} ***"
    folds = np.array([r["fold"] for r in rows], int)
    print("=" * 112)
    print(f"Q2  THE TORSION-SPACE LANDSCAPE OF THE TWO POTENTIALS   n = {n} targets{tag}")
    print("    BASIS: BUILT CHAIN.  Continuous torsion space, no lattice, no binary encoding.")
    print("=" * 112)
    print("\nGATE GC20b -- the fast chain-rule gradient against a finite-difference gradient")
    g = o.get("gate", {})
    print(f"    {json.dumps(g)[:600]}")

    print("\n1.  SCALE-INVARIANT LOCAL CURVATURE, at the pool-medoid start "
          "(a LOCAL statement -- see section 3 for the global one)")
    metrics = o["config"]["METRICS"]
    print(f"{'metric':<22s}{'Legacy':>12s}{'AMBER':>12s}{'paired diff':>34s}{'W/L':>10s}")
    for mname in metrics:
        a = np.array([r["start"]["legacy"][mname] for r in rows], float)
        b = np.array([r["start"]["amber"][mname] for r in rows], float)
        ok = np.isfinite(a) & np.isfinite(b)
        st = PP(b[ok], a[ok], folds=folds[ok])
        print(f"{mname:<22s}{np.nanmean(a):12.4f}{np.nanmean(b):12.4f}"
              f"{st['mean_diff']:+12.4f}[{st['ci95'][0]:+9.4f},{st['ci95'][1]:+9.4f}]"
              f"{st['n_better']:>5d}/{st['n_worse']:<4d}")
    print("    (AMBER minus Legacy.  Only scale-invariant metrics are compared across "
          "potentials -- F-C2b.)")
    print("\n    WITHIN-potential quantities, which HAVE UNITS and are never compared across:")
    for mname in ("grad_norm", "lam_max_abs", "E0"):
        a = np.array([r["start"]["legacy"].get(mname, np.nan) for r in rows], float)
        b = np.array([r["start"]["amber"].get(mname, np.nan) for r in rows], float)
        print(f"      {mname:<14s} Legacy {np.nanmean(a):14.4f}   AMBER {np.nanmean(b):14.4f}"
              f"   (kcal/mol-based units for AMBER, arbitrary weighted score for Legacy)")

    print("\n1b. IS EITHER POTENTIAL A FUNCTION OF THE STRUCTURE ALONE?  "
          "(phi[0] moves NO atom)")
    il = np.array([r["inert_max"]["legacy"] for r in rows], float)
    ia = np.array([r["inert_max"]["amber"] for r in rows], float)
    print(f"    |dE/dphi0|   AMBER  max over targets {ia.max():.3e}  "
          f"(EXACTLY zero on {int((ia == 0).sum())}/{len(ia)} targets)")
    print(f"                 Legacy max over targets {il.max():.3e}  "
          f"median {np.median(il):.3e}  (exactly zero on {int((il == 0).sum())}/{len(il)})")
    ef = np.array([r["legacy_channel"]["expl_frac"] for r in rows], float)
    cf = np.array([r["legacy_channel"]["coord_frac"] for r in rows], float)
    ad = np.array([r["legacy_channel"]["additivity_rel"] for r in rows], float)
    cc = np.array([r["legacy_channel"]["cos_coord_amber"] for r in rows], float)
    ct = np.array([r["legacy_channel"]["cos_total_amber"] for r in rows], float)
    print(f"    Legacy gradient channels: explicit-argument {ef.mean():.4f} of ||g|| "
          f"(median {np.median(ef):.4f}, max {ef.max():.4f}); coordinate {cf.mean():.4f}; "
          f"additivity residual max {ad.max():.2e}")
    gc = np.array([r["grad_cos"] for r in rows], float)
    st = PP(gc, np.zeros_like(gc), folds=folds)
    print(f"    cos(grad_Legacy, grad_AMBER) [SCALE-FREE]: {gc.mean():+.4f} "
          f"[{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}]  median {np.median(gc):+.4f}  "
          f"{int((gc > 0).sum())}/{len(gc)} positive")
    print(f"    cos(coordinate channel of Legacy, grad_AMBER): {np.nanmean(cc):+.4f}  "
          f"vs cos(total Legacy, grad_AMBER): {np.nanmean(ct):+.4f}  "
          f"(does dropping the structure-invisible channel bring Legacy closer to AMBER?)")

    print("\n2.  MULTI-START MINIMISATION IN TORSION SPACE  (identical starts, identical "
          "optimiser, identical budget; only H changes)")
    print(f"{'':<10s}{'E-drop':>10s}{'|g| drop':>10s}{'n distinct':>12s}"
          f"{'basin width':>13s}{'d(theta) moved':>15s}{'Ca-RMSD start':>15s}"
          f"{'Ca-RMSD end':>13s}")
    for pot in ("legacy", "amber"):
        f = lambda k: np.array([r["multi"][pot][k] for r in rows], float)   # noqa: E731
        print(f"{pot:<10s}{np.nanmean(f('e_drop_frac')):10.4f}"
              f"{np.nanmean(f('g_drop_frac')):10.4f}{np.nanmean(f('n_distinct')):12.2f}"
              f"{np.nanmean(f('basin_width')):13.4f}{np.nanmean(f('theta_moved')):15.4f}"
              f"{np.nanmean(f('rmsd_start')):15.4f}{np.nanmean(f('rmsd_end')):13.4f}")
    for pot in ("legacy", "amber"):
        s = np.array([r["multi"][pot]["rmsd_start"] for r in rows], float)
        e = np.array([r["multi"][pot]["rmsd_end"] for r in rows], float)
        print("    " + cmp_(e, s, folds, f"{pot}: end - start"))

    print("\n3.  BARRIERS AND CONNECTIVITY  (straight-line torus interpolation between the "
          "distinct minima found)")
    print("    READ THIS BEFORE THE NUMBERS.  A STRAIGHT LINE in torsion space between two")
    print("    peptide conformations passes through steric OVERLAP, so this is an UPPER BOUND on")
    print("    the true barrier, and for AMBER it is dominated by the r^-12 wall rather than by")
    print("    the landscape's connectivity.  It is reported as a MEDIAN (the mean is meaningless")
    print("    at these magnitudes) and it is a DIAGNOSTIC, not a barrier height: a real minimum")
    print("    energy path would go around, and this experiment does not compute one.")
    print(f"{'':<10s}{'n pairs':>9s}{'barrier/depth median':>22s}{'mean':>14s}"
          f"{'frac with barrier':>19s}")
    for pot in ("legacy", "amber"):
        f = lambda k: np.array([r["multi"][pot][k] for r in rows], float)   # noqa: E731
        br = f("barrier_rel")
        print(f"{pot:<10s}{np.nanmean(f('n_pairs')):9.2f}{np.nanmedian(br):22.4g}"
              f"{np.nanmean(br):14.4g}{np.nanmean(f('frac_barrier')):19.4f}")

    print("\n4.  *** THE TEST THAT DECIDES F-C2 ***  does any landscape metric PREDICT the "
          "final Ca-RMSD?")
    print(f"    POWER, STATED BEFORE THE TABLE: at n = {n} targets a Spearman needs "
          f"|rho| > {1.96/np.sqrt(max(n-3,1)):.2f} to clear zero at 95%.  A null here is "
          f"'NOT MEASURED for |rho| below that', NOT 'no relationship'.")
    print("    Spearman across targets between each metric and the Ca-RMSD the SAME optimiser")
    print("    reaches on the SAME potential from the SAME start, with a label-permutation null.")
    print("    TWO CLASSES, AND ONLY ONE OF THEM ANSWERS THE QUESTION (derive the operator before")
    print("    interpreting its statistic, BRIEF section 9):")
    print("      START-SIDE  (neg_frac, nearzero_frac, cond_med, gap_rel, part_ratio, aniso,")
    print("                   spec_skew, grad_norm) are properties of the landscape AT THE START,")
    print("                   computed before the optimiser runs.  These are PREDICTIVE.")
    print("      OUTCOME-SIDE (basin_width, n_distinct, barrier_rel) are computed FROM THE SAME")
    print("                   MINIMISATIONS whose RMSD is the endpoint -- 'the starts ended up far")
    print("                   apart' and 'the endpoint RMSD is high' are two readings of one fact.")
    print("                   These are NEARLY CIRCULAR and are labelled, never promoted.")
    print(f"{'metric':<24s}{'potential':>10s}{'class':>9s}{'rho':>8s}{'CI95':>22s}"
          f"{'perm null p':>13s}")
    hit = []
    for pot in ("legacy", "amber"):
        y = np.array([r["multi"][pot]["rmsd_end"] for r in rows], float)
        for mname in metrics + ["grad_norm", "basin_width", "n_distinct", "barrier_rel"]:
            if mname in metrics:
                x = np.array([r["start"][pot][mname] for r in rows], float)
            else:
                x = np.array([r["multi"][pot].get(mname, np.nan) for r in rows], float)
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 8 or np.nanstd(x[ok]) < 1e-12:
                continue
            rho, lo, hi, p = _rho_ci(x[ok], y[ok])
            side = "START" if mname in metrics or mname == "grad_norm" else "outcome"
            star = "  *" if (lo > 0 or hi < 0) else ""
            if star and side == "START":
                hit.append((pot, mname, round(rho, 3)))
            print(f"{mname:<24s}{pot:>10s}{side:>9s}{rho:8.3f}  [{lo:+8.3f},{hi:+8.3f}]"
                  f"{p:13.3f}{star}")
    print(f"\n    metrics with a CI excluding zero: {hit if hit else 'NONE'}")
    print(f"    F-C2 fires iff NO metric predicts Ca-RMSD  ->  "
          f"{'FIRES' if not hit else 'does NOT fire'}")

    fn_ = o.get("frame_null")
    if fn_:
        print(f"\nROTATED-FRAME NULL (exactly zero by construction; reported with its MAXIMUM): "
              f"{json.dumps(fn_)}")
    print("=" * 112)


def _rho_ci(x, y, n_boot=4000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(x)

    def sp(a, b):
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        ra -= ra.mean(); rb -= rb.mean()
        d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
        return float((ra * rb).sum() / d) if d > 0 else 0.0

    r0 = sp(x, y)
    bs = np.array([sp(x[i], y[i]) for i in (rng.integers(0, n, n) for _ in range(n_boot))])
    pn = np.array([sp(x, y[rng.permutation(n)]) for _ in range(2000)])
    p = float((np.abs(pn) >= abs(r0)).mean())
    return r0, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), p


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "q1"
    extra = sys.argv[2] if len(sys.argv) > 2 else None
    if which == "q1":
        q1()
    elif which == "pareto":
        pareto(extra or "c_pareto.json")
    elif which == "land":
        land(extra or "c_land.json")
    else:
        raise SystemExit("q1 | pareto | land")
