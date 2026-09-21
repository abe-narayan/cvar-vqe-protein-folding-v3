# -*- coding: utf-8 -*-
"""Recompute every headline number of Sprint 31 from its own artefact.

Built on `s30/s30_verify.py` (36 headline numbers, 0 mismatches), which this extends in the
four ways S31 asked for:

  1. every headline number is recomputed FROM ARTEFACTS, never read out of prose;
  2. EVERY PATH ANY LEDGER ENTRY CLAIMS TO HAVE WRITTEN must exist -- and here that is
     AUTOMATIC: the ledger is parsed for the paths it names, rather than a hand-kept list,
     because a hand-kept list is exactly what fails when someone forgets to add to it.
     S30 added this check after its own opening entry claimed to have saved the charter and
     had not -- the FIFTH instance in this project of prose naming a path that does not exist;
  3. where a number has TWO reporting bases (CA point cloud and built chain) both are checked,
     and any registered number whose basis is unstated is FLAGGED;
  4. S30's known sign-convention trap is carried forward and asserted:
     `s24.stats_lib.compare` is LOWER-IS-BETTER (d = a - b, negative = a better), so for a
     PREFERENCE RATE or any higher-is-better statistic its `verdict` string is INVERTED.
     QUOTE THE GATE, NEVER `.verdict`.

Usage:  python s31/s31_verify.py            # verify
        python s31/s31_verify.py --paths    # only the ledger-path audit
"""
import io
import json
import os
import re
import statistics as st
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OK, BAD, MISSING, FLAG = [], [], [], []

#: the two reporting bases this project's endpoint has.  A registered number must name one.
BASES = ("chain", "cloud", "selection", "in-band", "torsion", "none")


def load(path):
    if not os.path.exists(path):
        MISSING.append(path)
        return None
    try:
        return json.load(io.open(path, encoding="utf-8"))
    except Exception as e:
        MISSING.append("%s (%s)" % (path, e))
        return None


def dig(obj, *keys, **kw):
    cur = obj
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return kw.get("default")
        cur = cur[k]
    return cur


def check(label, claimed, actual, tol=5e-4, basis=None):
    """Recompute one headline number.  `basis` is MANDATORY for any RMSD-valued number."""
    if basis is None:
        FLAG.append(("basis unstated", label))
    elif basis not in BASES:
        FLAG.append(("basis not recognised: %s" % basis, label))
    if actual is None:
        MISSING.append(label)
        print("%-58s claim %-11s ARTEFACT KEY NOT FOUND" % (label[:58], claimed))
        return
    good = abs(float(claimed) - float(actual)) <= tol
    (OK if good else BAD).append((label, claimed, actual))
    print("%-58s claim %-11s actual %-11s %s  [%s]"
          % (label[:58], round(float(claimed), 5), round(float(actual), 5),
             "MATCH" if good else "*** MISMATCH ***", basis or "BASIS UNSTATED"))


def exact(label, claimed, actual):
    """A number that must match BIT-EXACTLY, not to a tolerance."""
    good = (actual is not None) and float(claimed) == float(actual)
    (OK if good else BAD).append((label, claimed, actual))
    print("%-58s claim %-11s actual %-11s %s"
          % (label[:58], claimed, actual, "EXACT" if good else "*** NOT EXACT ***"))


def show(label, value, note=""):
    print("%-58s %-24s %s" % (label[:58], str(value)[:24], note))


def rows(path, key, **flt):
    if not os.path.exists(path):
        MISSING.append(path)
        return None
    out = []
    for line in io.open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except Exception:
            continue
        if all(r.get(k) == v for k, v in flt.items()) and key in r:
            out.append(r[key])
    return out or None


def jrows(pattern, key, **flt):
    """Same, globbed over shard files."""
    import glob
    out = []
    for f in sorted(glob.glob(pattern)):
        got = rows(f, key, **flt)
        if got:
            out.extend(got)
    return out or None


print("=" * 96)
print("SPRINT 31 -- recomputing headline numbers from artefacts")
print("=" * 96)

paths_only = "--paths" in sys.argv

if not paths_only:
    # ================================================== lane D, D-B: the projection pin
    print()
    print("--- lane D, defect D-B (S31-L6): the projection, pinned ---")
    pin = load("s31/results/s31_D_projection_pin.json")
    if pin:
        show("n targets", dig(pin, "n"), "")
        # THE claim of the entry: bit-identical, not merely close
        exact("reprojecting the same cloud twice: max |dCA|", 0.0,
              dig(pin, "determinism", "max_coord_abs_diff_over_targets"))
        exact("reprojecting the same cloud twice: max |d rmsd|", 0.0,
              dig(pin, "determinism", "max_rmsd_abs_diff_over_targets"))
        b = dig(pin, "determinism", "all_bit_identical")
        (OK if b is True else BAD).append(("all 126 bit-identical", True, b))
        print("%-58s %s" % ("  all targets bit-identical",
                            "TRUE" if b else "*** FALSE ***"))
        # no RNG: asserted against the operator pin, which is the machine-readable claim
        rng = dig(pin, "operator_pin", "rng", default="")
        ok = "NONE" in str(rng)
        (OK if ok else BAD).append(("operator pin records rng=NONE", "NONE", rng))
        print("%-58s %s" % ("  operator pin: rng", rng))
        show("  operator pin digest", dig(pin, "operator_pin", "digest"),
             "changes if ANY projection constant changes")
        show("  n starts", dig(pin, "operator_pin", "n_starts"), "")

        # the endpoint reproduces BIT-FOR-BIT from the pinned inputs
        check("S31 reprojection vs s29 canonical: mean diff", 0.0,
              dig(pin, "vs_s29_prod_canonical", "mean_diff"), tol=1e-12, basis="chain")
        exact("  ... and max |diff| over 126 targets", 0.0,
              dig(pin, "vs_s29_prod_canonical", "max_abs_diff"))
        check("S31 reprojection mean (== canonical 3.2105)", 3.2105,
              dig(pin, "reprojection", "mean_chain_pass1"), basis="chain")
        check("  its CA cloud (the OTHER basis)", 3.0483,
              dig(pin, "reprojection", "mean_cloud"), basis="cloud")

        # the disagreement with S27, which is the defect itself
        check("vs s27 DIS: mean chain diff", -0.0021,
              dig(pin, "vs_s27_DIS", "mean_diff"), basis="chain")
        check("vs s27 DIS: max |chain diff| (2LNG)", 0.5174,
              dig(pin, "vs_s27_DIS", "max_abs_diff"), basis="chain")
        cd = dig(pin, "vs_s27_DIS", "cloud_max_abs_diff")
        ok = cd is not None and cd < 1e-12
        (OK if ok else BAD).append(("vs s27: the CLOUD is identical to 1e-12", "<1e-12", cd))
        print("%-58s %-11s %s" % ("  ... while the CLOUD agrees to", "%.2e" % (cd or 0),
                                  "MATCH" if ok else "*** MISMATCH ***"))
        show("  -> the cloud is exact; the CHAIN is what differs",
             "both bases checked", "this is the whole defect")

        # the mechanism
        print()
        print("--- the mechanism: conditioning, not stochasticity ---")
        fr = dig(pin, "flip_register") or {}
        show("lam=0 branch margin, median", "%.2e" % (fr.get("median_margin_lam0") or 0),
             "the argmin is decided at numerical noise")
        show("  targets with lam=0 margin < 1e-6",
             "%s / %s" % (fr.get("n_lam0_margin_below_1e-6"), dig(pin, "n")),
             "NOT one pathological target")
        show("  targets with EITHER margin < 1e-6",
             "%s / %s" % (fr.get("n_margin_below_1e-6"), dig(pin, "n")),
             "min(lam=0, final) -- a DIFFERENT count, not interchangeable")
        pb = dig(pin, "perturbation") or {}
        show("1e-14 relative input perturbation moves the input by",
             "%.2e A" % (pb.get("max_input_coord_shift_A") or 0), "")
        show("  ... and the OUTPUT chain by, at worst",
             "%.4f A (%s)" % (pb.get("max_abs_output_shift_A") or 0, pb.get("worst_target")),
             "amplification ~1e13")
        n0 = pb.get("n_output_unchanged_to_1e-9")
        ok = n0 == 0
        (OK if ok else BAD).append(("NOT ONE target unchanged to 1e-9", 0, n0))
        print("%-58s %-11s %s" % ("  targets whose output was unchanged to 1e-9", n0,
                                  "MATCH" if ok else "*** MISMATCH ***"))

        # the instrument spread -- the rule every lane is held to
        sp = dig(pin, "instrument_spread") or {}
        check("instrument spread across 5 records", 0.0107, sp.get("spread_A"),
              tol=5e-4, basis="chain")
        show("  -> ANY built-chain claim below this", "IS INSIDE THE NOISE",
             "confirmed effect for scale: 0.0221 A")

    # ================================================== lane D, D-B: recomputed from raw rows
    print()
    print("--- lane D: the five circulating 'production built chain' values, recomputed ---")
    for lbl, path, key, flt, claimed in (
            ("s29_O prod (CANONICAL)", "s29/results/s29_O_chain_rows.jsonl",
             "rmsd_chain", {"item": "prod"}, 3.2105),
            ("s27 chain_rows DIS", "s27/results/chain_rows.jsonl",
             "rmsd_chain", {"config": "DIS"}, 3.2126),
            ("s30_P PROD_chain", "s30/results/s30_P_chain_rows.jsonl",
             "PROD_chain", {}, 3.2126),
            ("s30_P rec_fit", "s30/results/s30_P_chain_rows.jsonl",
             "rec_fit", {}, 3.2041)):
        v = rows(path, key, **flt)
        if v:
            check("  %s" % lbl, claimed, st.mean(v), basis="chain")
            if len(v) != 126:
                FLAG.append(("n != 126 (%d)" % len(v), lbl))
    for lbl, path, key, flt, claimed in (
            ("s29_O prod CLOUD", "s29/results/s29_O_chain_rows.jsonl",
             "rmsd_cloud", {"item": "prod"}, 3.0483),
            ("s27 DIS CLOUD", "s27/results/chain_rows.jsonl",
             "rmsd_cloud", {"config": "DIS"}, 3.0483),
            ("s30_P PROD_cloud", "s30/results/s30_P_chain_rows.jsonl",
             "PROD_cloud", {}, 3.0483)):
        v = rows(path, key, **flt)
        if v:
            check("  %s" % lbl, claimed, st.mean(v), basis="cloud")
    show("  -> the CLOUD agrees across all three records", "3.0483",
         "while the CHAIN spans 3.2041-3.2148: BOTH BASES CHECKED")

    # ================================================== lane D: the branch experiment
    print()
    print("--- lane D, S31-L8: the branch experiment (PRE-REGISTERED) ---")
    br = load("s31/results/s31_D_branch.json")
    if br:
        pre = br.get("prereg")
        ok = bool(pre) and os.path.exists(pre)
        (OK if ok else BAD).append(("the prereg it names exists", pre, ok))
        print("%-58s %-24s %s" % ("  prereg", pre, "exists" if ok else "*** MISSING ***"))
        ic = br.get("identity_check_vs_s29_canonical") or {}
        exact("  PROD reproduces canonical, max |diff|", 0.0, ic.get("max_abs_diff"))
        (OK if ic.get("PASS") else BAD).append(("identity gate PASS", True, ic.get("PASS")))
        check("  PROD mean (must be the canonical endpoint)", 3.2105,
              dig(br, "means", "PROD"), basis="chain")
        check("  B4 mean", 3.2050, dig(br, "means", "B4"), basis="chain")
        check("  ORACLE_B4 mean  [ORACLE / NOT DEPLOYABLE]", 3.1168,
              dig(br, "means", "ORACLE_B4"), basis="chain")

        for key, claimed_eff, claimed_x, gate in (
                ("ORACLE_B4_vs_PROD", -0.0938, 2.92, "BETTER (>=1.0x MDE)"),
                ("B4_vs_PROD", -0.0055, 0.44, "NOT A RESULT (<0.7x MDE)")):
            c = br.get(key) or {}
            check("  %s effect" % key, claimed_eff, c.get("effect"), basis="chain")
            x = abs(c.get("effect", 0)) / c["mde"] if c.get("mde") else None
            check("  %s xMDE" % key, claimed_x, x, tol=5e-3, basis="chain")
            g = dig(br, "gates", key)
            okg = (g == gate)
            (OK if okg else BAD).append(("  %s GATE" % key, gate, g))
            print("%-58s %-24s %s" % ("    gate", g, "MATCH" if okg else "*** MISMATCH ***"))
            show("    SE / MDE / W-L-tied",
                 "%.4f / %.4f / %d-%d-%d" % (c.get("se", 0), c.get("mde", 0),
                                             c.get("n_better", 0), c.get("n_worse", 0),
                                             c.get("n_tied", 0)), "")

        # the standing rule: MDE = 2.8016 * SE, recomputed rather than trusted
        for key in ("ORACLE_B4_vs_PROD", "B4_vs_PROD"):
            c = br.get(key) or {}
            if c.get("se"):
                check("  %s MDE == 2.8016*SE" % key, 2.8016 * c["se"], c.get("mde"),
                      tol=1e-6, basis="none")

        # structural fact registered BEFORE the result
        od = br.get("objective_domination_check") or {}
        n = br.get("n")
        ok = od.get("n_B4_objective_le_PROD") == n
        (OK if ok else BAD).append(("B4 objective <= PROD on all n (BY CONSTRUCTION)",
                                    n, od.get("n_B4_objective_le_PROD")))
        print("%-58s %-24s %s" % ("  B4 objective <= PROD on all targets",
                                  "%s/%s" % (od.get("n_B4_objective_le_PROD"), n),
                                  "MATCH (structural, NOT evidence)" if ok else "*** MISMATCH ***"))

        cond = br.get("conditioning") or {}
        show("  conditioning: median margin lam=0 -> lam=0.3",
             "%.2e -> %.2e" % (cond.get("median_margin_PROD_lam0") or 0,
                               cond.get("median_margin_B4_lam") or 0),
             "%.0fx" % (cond.get("improvement_factor_median") or 0))
        show("  targets below 1e-6 margin: PROD -> B4",
             "%s -> %s" % (cond.get("n_PROD_margin_below_1e-6"),
                           cond.get("n_B4_margin_below_1e-6")), "")

        # the concentration rule: the median-vs-mean gap is NOT a verdict on its own
        cc = dig(br, "tail", "concentration") or {}
        show("  concentration: mean / median",
             "%.4f / %.4f" % (cc.get("mean") or 0, cc.get("median") or 0),
             "a gap is SUGGESTIVE, not a finding")
        show("  drop-top-10 vs its UNIFORM-EFFECT null",
             "%.4f at pct %.3f" % (cc.get("drop_top10_mean") or 0,
                                   cc.get("pctile_in_null") or 0),
             "flag=%s -> NOT ESTABLISHED" % cc.get("flag"))
        okf = (cc.get("flag") is False)
        (OK if okf else BAD).append(("concentration flag is False (null says ordinary)",
                                     False, cc.get("flag")))

    # ================================================== the sign-convention trap (from S30)
    print()
    print("--- the sign-convention trap, carried forward from S30 and re-asserted ---")
    print("  s24.stats_lib.compare is LOWER-IS-BETTER: d = a - b, negative = a better.")
    print("  For a PREFERENCE RATE or any higher-is-better statistic its .verdict is INVERTED.")
    print("  QUOTE THE GATE, NEVER .verdict.")
    try:
        sys.path.insert(0, ROOT)
        from s24 import stats_lib as ST
        import numpy as np
        # a is uniformly WORSE on a lower-is-better statistic
        # compare() REFUSES a verdict without fold-clustered CIs -- which is itself the
        # right behaviour and is asserted here, then folds are supplied for the real check.
        cnf = ST.compare(np.arange(10.0) + 1.0, np.arange(10.0), label="trap selftest nofold")
        okr = "NOT MEASURED" in str(cnf.get("verdict", ""))
        (OK if okr else BAD).append(("compare() refuses a verdict with no folds",
                                     "NOT MEASURED", cnf.get("verdict")))
        c = ST.compare(np.arange(10.0) + 1.0, np.arange(10.0),
                       folds=np.array([0, 1, 2, 3, 4] * 2), label="trap selftest")
        ok = c["effect"] > 0 and str(c.get("verdict", "")).upper().startswith("WORSE")
        (OK if ok else BAD).append(("compare() is lower-is-better (d=a-b)",
                                    "effect>0 & WORSE", "%.3f / %s"
                                    % (c["effect"], c.get("verdict"))))
        print("%-58s %-24s %s" % ("  selftest: a uniformly larger than b",
                                  "effect %+.3f, %s" % (c["effect"], c.get("verdict")),
                                  "CONVENTION HOLDS" if ok else "*** CONVENTION CHANGED ***"))
        check("  MDE_K is 2.8016", 2.8016, getattr(ST, "MDE_K", None), tol=1e-4, basis="none")
    except Exception as e:
        MISSING.append("s24.stats_lib selftest (%s)" % e)

    # ================================================== D-A and D-C, the code fixes
    print()
    print("--- lane D, S31-L7: the two code fixes, asserted against the SOURCE ---")
    try:
        src = io.open("core/pipeline.py", encoding="utf-8").read()
        bad = "that is the component's measured role" in src and "CORRECTED" not in src
        ok = ("S25-L5" in src and "-0.1405" in src and "0.68x MDE" in src
              and "NOT MEASURED" in src and not bad)
        (OK if ok else BAD).append(("D-A: the withdrawn +0.113 claim is corrected in place",
                                    "S25-L5 + replacement figure present", ok))
        print("%-58s %s" % ("  core/pipeline.py cites S25-L5 and -0.1405 at 0.68x MDE",
                            "FIXED" if ok else "*** STILL STALE ***"))
        okq = "Do not quote a number for the CVaR tail" in src
        (OK if okq else BAD).append(("D-A: the docstring forbids quoting a number", True, okq))
    except Exception as e:
        MISSING.append("core/pipeline.py (%s)" % e)

    try:
        jr = io.open("s26/jobrun.py", encoding="utf-8").read()
        gv = io.open("s26/governor.py", encoding="utf-8").read()
        derived = "_governor_band" in jr and "CPU_START = _GOV_CPU_CEILING" in jr
        (OK if derived else BAD).append(("D-C: jobrun derives its gate from governor.py",
                                         True, derived))
        print("%-58s %s" % ("  s26/jobrun.py imports the governor's constants",
                            "DERIVED" if derived else "*** STILL DUPLICATED ***"))
        # and they actually agree at runtime
        import importlib.util as ilu
        sp = ilu.spec_from_file_location("s31_gov_check", "s26/governor.py")
        gm = ilu.module_from_spec(sp)
        sp.loader.exec_module(gm)
        sp2 = ilu.spec_from_file_location("s31_jr_check", "s26/jobrun.py")
        jm = ilu.module_from_spec(sp2)
        sp2.loader.exec_module(jm)
        exact("  jobrun.CPU_START == governor.CPU_CEILING", gm.CPU_CEILING, jm.CPU_START)
        exact("  jobrun.CEILING == governor.CEILING - 1", gm.CEILING - 1.0, jm.CEILING)
        okc = jm._cap() <= gm.MAX_JOBS
        (OK if okc else BAD).append(("  jobrun cap clamped to governor.MAX_JOBS",
                                     "<= %d" % gm.MAX_JOBS, jm._cap()))
        show("  live cap / governor MAX_JOBS", "%s / %s" % (jm._cap(), gm.MAX_JOBS), "")
        # the launcher must NOT block in the band the charter asks for
        okband = not (95.0 > jm.CPU_START)
        (OK if okband else BAD).append(("  a launch is permitted at 94-95% CPU",
                                        "CPU_START >= 95", jm.CPU_START))
        print("%-58s %-24s %s" % ("  launch permitted in the charter's 94-95% CPU band",
                                  "CPU_START %.1f" % jm.CPU_START,
                                  "YES" if okband else "*** STILL BLOCKED ***"))
    except Exception as e:
        MISSING.append("s26 governor/jobrun agreement check (%s)" % e)

# ================================================== LANE F -- S31-L11 (F3) and the coh diagnostic
print()
print("--- lane F: S31-L11, the prefix-length order-statistic audit ---")
try:
    f3 = load("s31/results/s31_F3_prefix.json")
    check("F3a bestm128 vs production", -0.3079,
          dig(f3, "F3a_chain", "cmp", "effect"), basis="chain")
    check("F3a bestm128 mean", 2.9027, dig(f3, "F3a_chain", "cmp", "mean_a"), basis="chain")
    check("F3a production mean (re-projected by S29, same rows)", 3.2105,
          dig(f3, "F3a_chain", "cmp", "mean_b"), basis="chain")
    check("F3a bestm128 vs m=75", -0.2879, dig(f3, "F3a_cloud", "cmp", "effect"), basis="cloud")
    # the two registered bars, recomputed from the artefact rather than read from prose
    gp = dig(f3, "F3e_matched_random_family", "gain_prefix")
    gr = dig(f3, "F3e_matched_random_family", "gain_random_mean")
    check("F3e matched random family, gain", -0.4279, gr, basis="cloud")
    check("F3e share of the prefix gain", 1.4864, gr / gp, basis="cloud")
    exact("F3e BAR 'the m axis is an order statistic' FIRES (>= 0.80)",
          True, bool(gr / gp >= 0.80))
    exact("F3e prefix curve reproduces S29's curve128 bit-for-bit", 0.0,
          dig(f3, "F3e_matched_random_family", "prefix_curve_repro_maxdev"))
    check("F3c ORACLE global m", -0.0018,
          dig(f3, "F3c_global_m", "ORACLE_global", "cmp", "effect"), basis="cloud")
    check("F3c leave-fold-out m (WORSE)", 0.0079,
          dig(f3, "F3c_global_m", "LFO", "cmp", "effect"), basis="cloud")
    check("F3d split-half transfer, fraction of the oracle gain", 0.0123,
          dig(f3, "F3d_split_half_transfer", "frac_of_oracle"), basis="cloud")
    be = dig(f3, "F3f_summary", "best_effect")
    bm = dig(f3, "F3f_native_free_rules", dig(f3, "F3f_summary", "best_rule"), "cmp", "mde")
    check("F3f best native-free LFO m-rule", -0.0221, be, basis="cloud")
    exact("F3f BAR 'no part of it is deployable' FIRES (> -0.7x its own MDE)",
          True, bool(be > -0.7 * bm))
    # the order-statistic growth curve must be MONOTONE and UNSATURATED in both families
    for fam in ("prefix", "random"):
        g = dig(f3, "F3e_MECHANISM", "order_statistic_growth_gain_vs_k", fam)
        ks = sorted((int(k) for k in g), key=int)
        mono = all(g[str(ks[i + 1])] <= g[str(ks[i])] + 1e-12 for i in range(len(ks) - 1))
        still = (g[str(ks[-1])] - g[str(ks[-2])]) < -0.005
        exact("F3e growth curve MONOTONE in k (%s)" % fam, True, mono)
        exact("F3e growth curve STILL FALLING at k=128 (%s)" % fam, True, still)
except Exception as e:
    MISSING.append("lane F S31-L11 block (%s)" % e)

print()
print("--- lane F: the affine-hull / coh diagnostic (ORACLE) ---")
try:
    fc = load("s31/results/s31_F_coh.json")
    exact("coh(uniform mean in pair space) == 1 exactly (lane B's derivation)",
          True, abs(dig(fc, "coh", "uniform_mean_pairspace", "mean") - 1.0) < 1e-12)
    check("coh(coordinate average) reproduces lane B's 0.9780", 0.9780,
          dig(fc, "coh", "AVG", "mean"), basis="in-band")
    check("coh(ORACLE best member) reproduces lane B's 0.6708", 0.6708,
          dig(fc, "coh", "ORACLE_best", "mean"), basis="in-band")
    check("coh(AVG_SEP)", 0.9689, dig(fc, "coh", "AVG_SEP", "mean"), basis="in-band")
    check("coh(MED)", 0.8886, dig(fc, "coh", "MED", "mean"), basis="in-band")
    exact("AVG_SEP does NOT pass the 0.6931 admission bar",
          False, dig(fc, "coh", "AVG_SEP", "admitted"))
    exact("AVG_SEP passes on 0 of 126 targets", 0.0,
          dig(fc, "coh", "AVG_SEP", "frac_targets_under_bar"))
    check("AVG_SEP affine-hull residual, RMS per coordinate", 0.0452,
          dig(fc, "affine_hull_departure", "AVG_SEP_residual_rms_A", "mean"), basis="none")
except Exception as e:
    MISSING.append("lane F coh block (%s)" % e)

# ================================================== EVERY PATH THE LEDGER CLAIMS TO HAVE WRITTEN
print()
print("--- every path any S31 ledger entry names (parsed from the ledger, not hand-kept) ---")
PATH_RE = re.compile(r"`([A-Za-z0-9_./\\-]+\.(?:py|json|jsonl|md|npz|csv|txt|pptx|log))`")
BRACE_RE = re.compile(r"\{[0-9,]+\}")


def expand(tok):
    """`..._shard{0,1,2,3}of4.jsonl` -> the four real paths."""
    m = BRACE_RE.search(tok)
    if not m:
        return [tok]
    return [tok[:m.start()] + v + tok[m.end():] for v in m.group(0)[1:-1].split(",")]


claimed = []
for src in ("s31/LEDGER.md", "s31/MULTIPLICITY.md", "s31/STATE.md", "s31/S31_CONTRACT.md"):
    if not os.path.exists(src):
        MISSING.append(src)
        continue
    txt = io.open(src, encoding="utf-8").read()
    for tok in PATH_RE.findall(txt):
        for p in expand(tok):
            p = p.replace("\\", "/")
            #: only paths this sprint claims to have WRITTEN; a citation of an older sprint's
            #: artefact is checked too, because a dead citation is the same failure.
            if p.split("/")[0] in ("s31", "s30", "s29", "s27", "s26", "s25", "s24",
                                   "s20", "s12", "core", "tests", "verify"):
                claimed.append((src, p))

seen = {}
for src, p in claimed:
    seen.setdefault(p, set()).add(src)
nmiss = 0
for p in sorted(seen):
    good = os.path.exists(p)
    if good:
        OK.append(("path exists: %s" % p, "exists", "yes"))
    else:
        BAD.append(("path claimed but ABSENT: %s" % p, "exists", "*** MISSING ***"))
        nmiss += 1
        print("  *** MISSING *** %-56s  named in %s" % (p, ", ".join(sorted(seen[p]))))
print("  %d distinct paths named across the S31 documents; %d exist, %d MISSING"
      % (len(seen), len(seen) - nmiss, nmiss))
if nmiss == 0:
    print("  every path the ledger claims is on disk "
          "(the S30-L0 failure, which was the fifth of its kind, does not recur)")

# ================================================== basis audit
print()
print("--- basis audit: any registered number whose reporting basis is unstated ---")
if not FLAG:
    print("  none -- every checked number named one of %s" % (BASES,))
for why, lbl in FLAG:
    print("  FLAG %-28s %s" % (why, lbl))

print()
print("=" * 96)
print("MATCHED: %d    MISMATCHED: %d    KEYS/FILES NOT FOUND: %d    FLAGGED: %d"
      % (len(OK), len(BAD), len(MISSING), len(FLAG)))
for lbl, c, a in BAD:
    print("  MISMATCH  %s: recorded %s vs artefact %s" % (lbl, c, a))
for m in MISSING[:20]:
    print("  NOT FOUND %s" % m)
print("=" * 96)
sys.exit(1 if (BAD or MISSING) else 0)
