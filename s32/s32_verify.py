# -*- coding: utf-8 -*-
"""Recompute every headline number of Sprint 32 from the artefact that produced it.

Built on `s31/s31_verify.py` (271 checks), which is built on `s30/s30_verify.py`.  It keeps
all three of S31's guarantees and adds the two S32 earned in its first hours:

  1. every headline number is recomputed FROM ARTEFACTS, never read out of prose;
  2. every path any S32 document NAMES must exist -- parsed from the documents, not
     hand-kept, because a hand-kept list is exactly what fails when someone forgets it;
  3. CROSS-BASIS: a number is matched against BOTH reporting bases of every persisted ladder
     item and FLAGGED if it fits the basis it does NOT name.  (S31's worst defect: a built
     chain value 2.1435 labelled "CA cloud" in the headline, understating the sprint's own
     headroom by 0.16 A.  An audit that only asks whether a basis is NAMED cannot catch a
     MISSTATED one.)
  4. **NEW, S32-D1: the OBJECT audit.**  There are FIVE distinct chain-like means over the
     same 126 targets, spanning 0.031 A, and the canonical endpoint is NONE of the four
     cached scalars.  Any number attributed to the wrong object is flagged even when its
     basis word ("chain") is correct.  A basis is not an object.
  5. **NEW, S32-D1: the BIT-IDENTITY rule.**  The projection is deterministic and
     DISCONTINUOUS: 7e-15 A of input change gives 0.10-0.15 A of output change.  Any chain
     contrast quoted across jobs must carry the bit-identity check; this scans the documents
     for cross-job chain claims that do not.

SELF-TESTS.  Contract rule 5: a verification must be able to fail.  Each of the three audits
above ships with a self-test that feeds it the real historical defect and asserts it trips,
AND feeds it the corrected version and asserts it does not.  If a self-test fails, this
script exits non-zero even when every number matches.

Usage:  python s32/s32_verify.py              # everything
        python s32/s32_verify.py --paths      # only the document-path audit
        python s32/s32_verify.py --selftest   # only the self-tests
"""
from __future__ import annotations

import glob
import io
import json
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OK, BAD, MISSING, FLAG = [], [], [], []
EXEMPT = []          # (file, line, audit, why) -- struck claims whose replacement IS stated

#: what each phrase audit accepts as "the replacement is stated in the same block".  The
#: exemption is a POSITIVE OBLIGATION, not a loophole: a strike with no replacement still flags.
A5_CORRECTION = ["total suppression", "gain_out", "0 on its", "cleans up", "hull floor",
                 "1.8290", "exactly 0"]
A6_CORRECTION = ["5.25", "5.2", "rank(aff", "32.88", "globally"]
A7_CORRECTION = ["0.3", "lam=0.3"]

#: the reporting bases this project's numbers have.  A registered RMSD number must name one.
BASES = ("chain", "cloud", "set-mean", "member", "selection", "in-band", "torsion", "none")

#: S32 documents scanned for paths, numbers and unsupported cross-job chain claims.
DOCS = ("s32/LEDGER.md", "s32/STATE.md", "s32/MULTIPLICITY.md", "s32/S32_CONTRACT.md",
        "s32/REPORT_S32.md", "s32/CAUSAL_MAP.md", "s32/THEORY_Q.md",
        "s32/PREREG_S32_D.md", "s32/PREREG_S32_P.md", "s32/PREREG_S32_Q.md",
        "s32/PREREG_S32_R.md", "s32/AUDIT_V.md")

PROD_KEY = "1fc9f2dcf489e2fb"
PROD_CACHE = os.path.join("bench_results", "cache", PROD_KEY)


# --------------------------------------------------------------------------- primitives
def load(path):
    if not os.path.exists(path):
        MISSING.append(path)
        return None
    try:
        return json.load(io.open(path, encoding="utf-8"))
    except Exception as e:                                      # noqa: BLE001
        MISSING.append("%s (%s)" % (path, e))
        return None


def dig(obj, *keys, **kw):
    cur = obj
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return kw.get("default")
        cur = cur[k]
    return cur


def jsonl(pattern, **flt):
    """Every row of every file matching `pattern`, filtered, DEDUPED by (pdb, item)."""
    out = {}
    for f in sorted(glob.glob(pattern)):
        for line in io.open(f, encoding="utf-8"):
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except Exception:                                   # noqa: BLE001
                continue
            if all(r.get(k) == v for k, v in flt.items()):
                out[(r.get("pdb"), r.get("item"))] = r
    return list(out.values())


def mean(xs):
    xs = [float(x) for x in xs if x is not None and x == x]
    return sum(xs) / len(xs) if xs else None


CHECKS = []          # (label, claimed, basis, object_name) -- consumed by the two audits


def check(label, claimed, actual, tol=5e-4, basis=None, obj=None):
    """Recompute one headline number.  `basis` is MANDATORY for any RMSD-valued number;
    `obj` names WHICH of the five chain-like objects it is, when basis == 'chain'."""
    CHECKS.append((label, float(claimed), basis, obj))
    if basis is None:
        FLAG.append(("basis unstated", label))
    elif basis not in BASES:
        FLAG.append(("basis not recognised: %s" % basis, label))
    elif basis == "chain" and obj is None:
        FLAG.append(("chain number with no OBJECT named (S32-D1)", label))
    if actual is None:
        MISSING.append(label)
        print("%-56s claim %-11s ARTEFACT KEY NOT FOUND" % (label[:56], claimed))
        return
    good = abs(float(claimed) - float(actual)) <= tol
    (OK if good else BAD).append((label, claimed, actual))
    print("%-56s claim %-11s actual %-11s %s  [%s]"
          % (label[:56], round(float(claimed), 6), round(float(actual), 6),
             "MATCH" if good else "*** MISMATCH ***", basis or "BASIS UNSTATED"))


def exact(label, claimed, actual):
    good = (actual is not None) and float(claimed) == float(actual)
    (OK if good else BAD).append((label, claimed, actual))
    print("%-56s claim %-11s actual %-11s %s"
          % (label[:56], claimed, actual, "EXACT" if good else "*** NOT EXACT ***"))


def show(label, value, note=""):
    print("%-56s %-26s %s" % (label[:56], str(value)[:26], note))


# ============================================ RETRACTION CONTEXTS (contract rule 13 vs the audits)
#: Rule 13 retracts IN PLACE and never deletes, so a phrase-matching audit fires forever on every
#: correctly-handled retraction, and the permanent-flag count grows monotonically until a genuine
#: new flag is invisible.  The exemption below is deliberately NOT "the word 'corrected' is near
#: it" -- that would let any claim be laundered by writing 'corrected' beside it.  A struck claim
#: is exempt ONLY IF ITS REPLACEMENT IS STATED IN THE SAME BLOCK.  A strike with no replacement is
#: itself a defect this project has committed, and is now caught rather than excused.
RETRACT_OPEN = re.compile(r"(?i)(~~|\bSTRUCK\b|\bRETRACT(ED|ION)?\b|\[corrected|\bWITHDRAWN\b|"
                          r"original wording|left standing per rule 13|\[false\b|is FALSE\b|"
                          r"\bbackwards\b|\bI wrote\b|first wrote|\bthe correction\b|"
                          r"\brefuted\b|\bwas wrong\b|\bannotated in place\b)")


def _blocks(raw):
    """Contiguous non-blank line runs, as (start, end_exclusive)."""
    out, s = [], None
    for k, l in enumerate(raw):
        if l.strip() and s is None:
            s = k
        elif not l.strip() and s is not None:
            out.append((s, k)); s = None
    if s is not None:
        out.append((s, len(raw)))
    return out


def retraction_exempt(raw, ln0, correction_tokens):
    """True iff line `ln0` sits in a retraction context AND that context states the replacement.

    `correction_tokens`: any one of these appearing in the same block (or inside the same ~~...~~
    span) discharges the retraction.  Empty list => a retraction context alone never exempts.
    """
    line = raw[ln0]
    blk = next(((a, b) for a, b in _blocks(raw) if a <= ln0 < b), (ln0, ln0 + 1))
    ctx = "\n".join(raw[blk[0]:blk[1]])
    struck = bool(RETRACT_OPEN.search(ctx)) or bool(re.search(r"~~.*~~", line))
    if not struck:
        return False, "not a retraction context"
    if not correction_tokens:
        return False, "retraction context, but no replacement token configured"
    low = ctx.lower()
    for t in correction_tokens:
        if t.lower() in low:
            return True, "struck, replacement '%s' stated in the same block" % t
    return False, "STRUCK WITH NO REPLACEMENT STATED IN THE BLOCK"


paths_only = "--paths" in sys.argv
selftest_only = "--selftest" in sys.argv

print("=" * 98)
print("SPRINT 32 -- recomputing headline numbers from artefacts")
print("=" * 98)


# ============================================================ THE FIVE OBJECTS (S32-L4/D1)
#: Recomputed here from the production cache and the s29 rows, so the L4 table cannot rot.
def prod_cache_means():
    keys = ("rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full", "top_m_mean",
            "pool_best", "top_m_best", "shipped", "fold")
    acc = {k: [] for k in keys}
    fs = sorted(glob.glob(os.path.join(PROD_CACHE, "*.json")))
    for f in fs:
        d = json.load(io.open(f, encoding="utf-8"))
        for k in keys:
            acc[k].append(d.get(k))
    return {k: mean(v) for k, v in acc.items() if k != "fold"}, len(fs), acc["fold"]


CACHE, N_CACHE, FOLDVEC = prod_cache_means()
LADDER_PROD = jsonl("s29/results/s29_O_chain_rows*.jsonl", item="prod")

#: name -> (value, source, one-line definition).  THE object registry.
OBJECTS = {
    "cloud (avg_ca)": (CACHE.get("rmsd_avg"), "bench_results/cache/%s/*.json:rmsd_avg" % PROD_KEY,
                       "top-75 coordinate average, unprojected"),
    "chain lam=0 (fit_ca)": (CACHE.get("rmsd_fit"), "bench_results/cache/%s/*.json:rmsd_fit" % PROD_KEY,
                             "projection arm with the Ramachandran penalty OFF"),
    "chain lam=0.3 EMITTED (ca)": (CACHE.get("rmsd_arm"), "bench_results/cache/%s/*.json:rmsd_arm" % PROD_KEY,
                                   "the chain the production pipeline itself emits"),
    "chain lam=0.3 + AMBER (amber_ca)": (CACHE.get("rmsd_full"), "bench_results/cache/%s/*.json:rmsd_full" % PROD_KEY,
                                         "the emitted chain after AMBER relaxation"),
    "chain CANONICAL (re-projection)": (mean([r["rmsd_chain"] for r in LADDER_PROD]),
                                        "s29/results/s29_O_chain_rows*.jsonl item=prod",
                                        "s12.instrument.project re-applied to s29_O_structs['prod']"),
    "set mean (top-75 members)": (CACHE.get("top_m_mean"), "bench_results/cache/%s/*.json:top_m_mean" % PROD_KEY,
                                  "mean ORACLE rr over the 75 selected members"),
}

if not selftest_only:
    print()
    print("--- S32-L4: the five chain-like objects, each recomputed from its own source ---")
    for nm, (v, src, d) in OBJECTS.items():
        print("  %-34s %s   %-52s  %s" % (nm, ("%.6f" % v) if v is not None else "  MISSING   ", src, d))
    _vals = [v for v in (OBJECTS[k][0] for k in OBJECTS) if v is not None]
    show("  spread of the four chain-like means",
         "%.4f A" % (max(v for k, v in ((k, OBJECTS[k][0]) for k in OBJECTS) if "chain" in k)
                     - min(v for k, v in ((k, OBJECTS[k][0]) for k in OBJECTS) if "chain" in k)),
         "none of these is a substitute for another")
    check("n targets in the production cache", 126, N_CACHE, tol=0, basis="none")
    check("n canonical prod chain rows", 126, len(LADDER_PROD), tol=0, basis="none")

    print()
    print("--- the endpoint, the two diagnostics, and the projection price ---")
    check("THE ENDPOINT (canonical re-projection)", 3.2105,
          OBJECTS["chain CANONICAL (re-projection)"][0], basis="chain",
          obj="chain CANONICAL (re-projection)")
    check("CA point cloud (diagnostic)", 3.0483, OBJECTS["cloud (avg_ca)"][0], basis="cloud")
    check("set mean over the top-75 (a third object)", 3.5507,
          OBJECTS["set mean (top-75 members)"][0], basis="set-mean")
    check("the chain PRODUCTION ITSELF EMITS (not the endpoint)", 3.2148,
          OBJECTS["chain lam=0.3 EMITTED (ca)"][0], basis="chain",
          obj="chain lam=0.3 EMITTED (ca)")
    check("endpoint minus emitted chain", -0.0043,
          (OBJECTS["chain CANONICAL (re-projection)"][0]
           - OBJECTS["chain lam=0.3 EMITTED (ca)"][0]), basis="chain",
          obj="chain CANONICAL (re-projection)")
    check("projection price, PRODUCTION (chain - cloud)", 0.1622,
          (OBJECTS["chain CANONICAL (re-projection)"][0] - OBJECTS["cloud (avg_ca)"][0]),
          basis="chain", obj="chain CANONICAL (re-projection)")

    # ---------------------------------------------------- lane V, step 4: the rebuild
    print()
    print("--- lane V, charter Step 4: the independent rebuild (s32_V_step4_endpoint.json) ---")
    s4 = load("s32/results/s32_V_step4_endpoint.json")
    if s4:
        by = {c["name"]: c for c in dig(s4, "checks", default=[]) if isinstance(c, dict)}
        check("rebuild: n targets", 126, dig(s4, "n"), tol=0, basis="none")
        check("rebuild: CA point cloud", 3.0483, dig(by, "CA point cloud mean (diagnostic)", "got"),
              basis="cloud")
        check("rebuild: set mean", 3.5507, dig(by, "set mean over top-75 members", "got"), basis="set-mean")
        check("rebuild: pool best K=500 [ORACLE]", 1.7108,
              dig(by, "pool best (K=500) [ORACLE]", "got"), tol=2e-3, basis="member")
        check("rebuild: top-75 best [ORACLE]", 2.3062,
              dig(by, "top-75 best [ORACLE]", "got"), tol=2e-3, basis="member")
        check("rebuild: shipped argmin [ORACLE]", 3.4540,
              dig(by, "shipped argmin [ORACLE label]", "got"), tol=2e-3, basis="member")
        exact("rebuild: top-75 set reproduced on N/126", 126, dig(s4, "n_sub_set_match"))
        v = dig(s4, "max_dev_avg_vs_stored")
        ok = v is not None and v < 1e-12
        (OK if ok else BAD).append(("coordinate average reproduces to <1e-12", "<1e-12", v))
        print("%-56s %-11s %s" % ("  |C_rebuilt - stored avg_ca| max", "%.2e" % (v or 0),
                                  "MATCH" if ok else "*** MISMATCH ***"))
        fc = dig(s4, "fold_counts") or {}
        ok = sorted(int(x) for x in fc.values()) == [23, 23, 25, 25, 30]
        (OK if ok else BAD).append(("fold sizes 25/23/25/23/30", [25, 23, 25, 23, 30], sorted(fc.values())))
        print("%-56s %-26s %s" % ("  fold structure", fc, "MATCH" if ok else "*** MISMATCH ***"))

    # ------------------------------------------- lane V: is the endpoint bit-reproducible?
    print()
    print("--- lane V, S32-D1: the endpoint's bit-reproducibility (s32_V_chain_bitexact.json) ---")
    bx = load("s32/results/s32_V_chain_bitexact.json")
    if bx:
        check("canonical mean recomputed from the s29 rows", 3.2105,
              dig(bx, "canonical_mean_s29"), basis="chain", obj="chain CANONICAL (re-projection)")
        a1 = dig(bx, "arms", "s29 struct npz ['prod']") or {}
        exact("re-projecting s29_O_structs['prod']: bit-identical on N/126",
              a1.get("n"), a1.get("n_bit_identical"))
        for nm in ("production cache avg_ca", "recomputed coordinate average"):
            a = dig(bx, "arms", nm) or {}
            if a:
                show("  SAME average, other float64 bits: %s" % nm,
                     "mean %+.6f  max|d| %.4f" % (a.get("mean_diff", 0), a.get("max_abs_diff", 0)),
                     "bit-identical on %s/%s" % (a.get("n_bit_identical"), a.get("n")))
        show("  max |cache avg_ca - s29 struct prod|",
             "%.2e A" % (dig(bx, "max_abs_cloud_diff_over_targets") or 0),
             "-> up to 0.52 A of built chain")
        st = dig(bx, "selftest_detects_ulp_perturbed_input")
        (OK if st else BAD).append(("bitexact selftest: a ULP-perturbed input IS detected", True, st))
        print("%-56s %s" % ("  selftest: ULP-perturbed input detected", "YES" if st else "*** NO ***"))

    # --------------------------------------- lane V: the ladder's narrowing order statistic
    print()
    print("--- lane V, S32-D2: the ladder's narrowing increments (s32_V_ladder_orderstat.json) ---")
    os_ = load("s32/results/s32_V_ladder_orderstat.json")
    if os_:
        inc = dig(os_, "increments") or {}
        check("K=500 -> 128 TOTAL [ORACLE, member basis]", 0.4350,
              inc.get("K=500 -> 128 TOTAL"), basis="member")
        check("  ... of which set-size order statistic", 0.2477,
              inc.get("  size effect (random 128)"), basis="member")
        check("  ... of which the DIS ordering (WRONG SIGN)", 0.1872,
              inc.get("  ordering effect (DIS 128 - random 128)"), basis="member")
        check("128 -> 75 TOTAL [ORACLE, member basis]", 0.1604,
              inc.get("128 -> 75 TOTAL"), basis="member")
        check("  ... of which set-size order statistic", 0.0907,
              inc.get("  size effect (random 75 of 128)"), basis="member")
        check("  ... of which the DIS ordering (WRONG SIGN)", 0.0696,
              inc.get("  ordering effect (DIS 75 - random 75)"), basis="member")
        show("  draws per target", dig(os_, "n_draws"), "contract rule 10: a draw control needs its own distribution")
        show("  recall of the oracle-best member into the top-128",
             "%.3f" % (dig(os_, "recall", "best500_in_top128") or 0),
             "vs 128/500 = 0.256 at random -- the score DOES find it more often")
    sp = load("s32/results/s32_V_top128_spread.json")
    if sp:
        pm = mean([v["pool_mean"] for v in sp.values()])
        tm = mean([v["top128_mean"] for v in sp.values()])
        ps = mean([v["pool_sd"] for v in sp.values()])
        ts = mean([v["top128_sd"] for v in sp.values()])
        p5 = mean([v["pool_p5"] for v in sp.values()])
        t5 = mean([v["top128_p5"] for v in sp.values()])
        check("mechanism: pool mean rr", 4.4533, pm, tol=2e-3, basis="member")
        check("mechanism: DIS top-128 mean rr", 3.5847, tm, tol=2e-3, basis="member")
        check("mechanism: pool sd rr", 1.3206, ps, tol=2e-3, basis="member")
        check("mechanism: DIS top-128 sd rr", 0.6384, ts, tol=2e-3, basis="member")
        check("mechanism: pool 5th pct rr", 2.6098, p5, tol=2e-3, basis="member")
        check("mechanism: DIS top-128 5th pct rr (UNCHANGED)", 2.6184, t5, tol=2e-3, basis="member")


    # ------------------------------------------- lane V: the ULP distribution of the endpoint
    print()
    print("--- lane V, S32-D1: the arithmetic-noise distribution of the endpoint ---")
    up = load("s32/results/s32_V_ulp_distribution.json")
    if up:
        check("ULP job: n targets", 126, dig(up, "n_targets"), tol=0, basis="none")
        check("ULP job: complete draws", 5, dig(up, "draws_complete"), tol=0, basis="none")
        v = dig(up, "max_abs_cloud_rmsd_shift")
        ok = v is not None and v < 1e-9
        (OK if ok else BAD).append(("the perturbation leaves the CLOUD unchanged", "<1e-9", v))
        print("%-56s %-11s %s" % ("  perturbation's max cloud-RMSD shift", "%.2e" % (v or 0),
                                  "MATCH" if ok else "*** MISMATCH ***"))
        check("draw-to-draw sd of the endpoint (THE RULE-20 NUMBER)", 0.0030,
              dig(up, "draw_sd"), tol=5e-4, basis="chain", obj="chain CANONICAL (re-projection)")
        show("  per-target |chain - canonical|",
             "mean %.4f p90 %.4f max %.4f" % (dig(up, "per_target_abs_delta", "mean") or 0,
                                              dig(up, "per_target_abs_delta", "p90") or 0,
                                              dig(up, "per_target_abs_delta", "max") or 0),
             "rule 3's floor: 0.0134 / 0.0329 / 0.2285")

    # ------------------------------------ lane V: is lane R's cos_align independent evidence?
    print()
    print("--- lane V: lane R's cos_align -- identity check (s32_V_cos_identity.json) ---")
    ci = load("s32/results/s32_V_cos_identity.json")
    if ci:
        v = dig(ci, "price_reconstructed_from_cos_max_abs_err")
        exact("  price rebuilt from (e,d,cos) alone: max abs error", 0.0, v)
        show("  -> cos_align is a bijection with the price given (e,d)",
             "NOT independent evidence", "quoting both double-counts one measurement")
        show("  cos_algebraic vs cos_direct (common frame)",
             "mean diff %.4f, max %.4f" % (dig(ci, "cos_difference", "abs_mean") or 0,
                                           dig(ci, "cos_difference", "max") or 0),
             "the Kabsch triangle IS near-Euclidean here")
        show("  production cos_align, mean", "%+.4f" % (dig(ci, "cos_algebraic", "mean") or 0),
             "registered null is 0 (orthogonal); negative = slightly anti-aligned")

    # ------------------------------- lane V: the adversarial replication of lane D's D1-T
    print()
    print("--- lane V: lane D's D1-T sign transfer, replicated with the controls it lacks ---")
    ds = load("s32/results/s32_V_D_signadversary.json")
    if ds:
        check("replication of lane D's DIS transfer (+0.1890)", 0.1890,
              dig(ds, "scorers", "DIS", "transfer"), tol=0.01, basis="in-band")
        for nm in ("DIS", "TYPICALITY", "RG", "NOISE", "DIS_DEDUP", "RG_DEDUP"):
            r = dig(ds, "scorers", nm)
            if r:
                print("  %-18s transfer %+0.4f  nullPERM %+0.4f  nullXTGT %+0.4f  globalSGN %+0.4f  %s"
                      % (nm, r["transfer"], r["null_perm"], r["null_xtarget_mean"],
                         r["global_sign_lfo"], r["verdict_per_target"]))
        nz = dig(ds, "scorers", "NOISE", "verdict_per_target")
        ok = nz == "NOT A RESULT"
        (OK if ok else BAD).append(("the pure-noise falsifier reports NOT A RESULT", "NOT A RESULT", nz))
        print("%-56s %s" % ("  falsifier: pure noise must be NOT A RESULT", nz))
        show("  band duplicates (the split-half independence assumption)",
             "%.1f%% of 75, on %d/126" % (100 * (dig(ds, "band_duplicates", "mean_frac") or 0),
                                          dig(ds, "band_duplicates", "n_targets_with_any") or 0),
             "dedup moves the transfer by <=0.012")

    # ---------------------------------------------------------------- lane Q: the readout identity
    print()
    print("--- lane Q, S32-L5: the readout program's gain (s32_Q1_sufficiency.json) ---")
    q1 = load("s32/results/s32_Q1_sufficiency.json")
    if q1:
        show("  label", dig(q1, "label"), "")
        show("  basis", dig(q1, "basis"), "")
        check("Q1-T1 identity, max relative error", 0.0, dig(q1, "T1", "err_identity_rel_max"),
              tol=1e-9, basis="none")
        check("Q1-T2 gain INSIDE the active affine hull", 1.0, dig(q1, "T2", "gain_in_mean"),
              tol=1e-4, basis="none")
        check("Q1-T2 gain OUTSIDE it (total suppression, not none)", 0.0,
              dig(q1, "T2", "gain_out_mean"), tol=1e-4, basis="none")
        check("active affine-hull dimension |S|-1, mean", 5.2540,
              dig(q1, "T2", "support_aff_dim_mean"), tol=1e-3, basis="none")
        check("ambient dimension d = 3n, mean", 38.8810, dig(q1, "T1", "d_mean"), tol=1e-3, basis="none")
        check("hull floor d (kabsch) [ORACLE]", 1.8290, dig(q1, "hull", "kabsch_mean"), basis="cloud")
        nz = dig(q1, "noise_ORACLE") or []
        for r in nz:
            dd = r["proj_mean"] - r["direct_mean"]
            print("  %-30s direct %.4f  projected %.4f   %s"
                  % ("noise eps=%.1f [ORACLE]" % r["eps"], r["direct_mean"], r["proj_mean"],
                     "the projection CLEANS UP %.4f A" % -dd if dd < 0
                     else "the projection costs %.4f A" % dd))

        # --- the contradiction audit: the phrase vs the artefact that is supposed to support it
        print()
        print("--- AUDIT 5: 'no noise suppression' vs the measured gain (S32-V, lane Q) ---")
        gout = dig(q1, "T2", "gain_out_mean")
        dim_out = (dig(q1, "T1", "d_mean") or 0) - (dig(q1, "T2", "support_aff_dim_mean") or 0)
        cleaned = [r for r in nz if r["proj_mean"] < r["direct_mean"] - 1e-6]
        contradicted = (gout is not None and gout < 1e-3 and dim_out > 1) or bool(cleaned)
        hits = []
        for src in DOCS:
            if not os.path.exists(src):
                continue
            raw = io.open(src, encoding="utf-8").read().splitlines()
            for ln0, line in enumerate(raw):
                if not (re.search(r"no\s+noise[- ]suppression", line, re.I) and contradicted):
                    continue
                ex, why = retraction_exempt(raw, ln0, A5_CORRECTION)
                if ex:
                    EXEMPT.append((src, ln0 + 1, "AUDIT 5", why))
                    continue
                hits.append((src, ln0 + 1, line.strip()[:80]
                             + ("   [%s]" % why if "NO REPLACEMENT" in why else "")))
        for s, ln, t in hits:
            FLAG.append(("'no noise suppression' but gain_out=%.1e on a %.0f-dim complement, "
                         "and %d noise rows ARE cleaned up" % (gout or 0, dim_out, len(cleaned)),
                         "%s:%d" % (s, ln)))
            print("  FLAG %s:%d  %s" % (s, ln, t))
        print("  gain_out_mean %.2e on a %.1f-dimensional complement; %d of %d noise rows show the "
              "projection CLEANING UP the estimate; %d unqualified claim(s) of 'no noise suppression'"
              % (gout or 0, dim_out, len(cleaned), len(nz), len(hits)))

        # --- the two-quantities-one-symbol audit: |S|-1 is 5.25, rank(aff{W}) is 32.88
        print()
        print("--- AUDIT 6: |S|-1 (the ACTIVE support, 5.25) vs rank(aff{W}) (ALL candidates, 32.88) ---")
        sdim = dig(q1, "T2", "support_aff_dim_mean")
        radim = dig(q1, "T1", "rank_aff_mean")
        n6 = 0
        for src in DOCS:
            if not os.path.exists(src):
                continue
            raw = io.open(src, encoding="utf-8").read().splitlines()
            for ln0, line in enumerate(raw):
                if not (re.search(r"\|S\|\s*[-−]\s*1", line) and re.search(r"\b3[0-9]\b", line)
                        and "rank" not in line.lower() and "globally" not in line.lower()):
                    continue
                ex, why = retraction_exempt(raw, ln0, A6_CORRECTION)
                if ex:
                    EXEMPT.append((src, ln0 + 1, "AUDIT 6", why))
                    continue
                FLAG.append(("'|S|-1' given as ~3x while |S|-1 = %.2f and rank(aff{W}) = %.2f"
                             % (sdim or 0, radim or 0), "%s:%d" % (src, ln0 + 1)))
                print("  FLAG %s:%d  %s %s" % (src, ln0 + 1, line.strip()[:78],
                                               "[%s]" % why if "NO REPLACEMENT" in why else ""))
                n6 += 1
        print("  |S|-1 mean %.2f (the ACTIVE support); rank(aff{W}) mean %.2f (ALL 128 candidates); "
              "%d conflation(s)" % (sdim or 0, radim or 0, n6))

    # --- AUDIT 7: the deployed projection arm is lambda=0.3, asserted from the SOURCE
    print()
    print("--- AUDIT 7: which lambda arm does the deployed chain come from? (asserted from source) ---")
    _src = io.open(os.path.join(ROOT, "s12", "instrument.py"), encoding="utf-8").read()
    _m = re.search(r"def project\(C, seq, fold, lam=([0-9.]+)", _src)
    _lam = float(_m.group(1)) if _m else None
    _returns_arm = bool(re.search(r'path\[0\.0\].*?path\[lam\]', _src, re.S)) and '"ca": np.asarray(arm[0]' in _src
    show("  s12/instrument.project default lam", _lam, "the deployed chain is path[lam], not path[0.0]")
    (OK if (_lam == 0.3 and _returns_arm) else BAD).append(
        ("instrument.project emits the lam=0.3 arm as 'ca'", "0.3 / arm", "%s / %s" % (_lam, _returns_arm)))
    n7 = 0
    for src in DOCS:
        if not os.path.exists(src):
            continue
        raw = io.open(src, encoding="utf-8").read().splitlines()
        for ln0, line in enumerate(raw):
            if not re.search(r"(final rung|last rung|the chain|deployed).{0,60}(runs at|is at|uses)\s*"
                             r"\*{0,2}\s*(λ|lam(bda)?)\s*=\s*0(?![.\d])", line, re.I):
                continue
            ex, why = retraction_exempt(raw, ln0, A7_CORRECTION)
            if ex:
                EXEMPT.append((src, ln0 + 1, "AUDIT 7", why))
                continue
            FLAG.append(("says the final rung runs at lambda=0; the source default is %s" % _lam,
                         "%s:%d" % (src, ln0 + 1)))
            print("  FLAG %s:%d  %s %s" % (src, ln0 + 1, line.strip()[:78],
                                           "[%s]" % why if "NO REPLACEMENT" in why else ""))
            n7 += 1
    print("  %d document line(s) still asserting lambda=0 for the deployed rung" % n7)


# ================================================ EVERY PATH ANY S32 DOCUMENT NAMES
if not selftest_only:
    print()
    print("--- every path any S32 document names (parsed from the documents, not hand-kept) ---")
PATH_RE = re.compile(r"`([A-Za-z0-9_./\\{},*-]+\.(?:py|json|jsonl|md|npz|csv|txt|pptx|log))`")
BRACE_RE = re.compile(r"\{[0-9,]+\}")
SPRINT_DIRS = ("s32", "s31", "s30", "s29", "s28", "s27", "s26", "s25", "s24", "s20", "s12",
               "s8", "core", "tests", "verify", "bench_results")


def expand(tok):
    m = BRACE_RE.search(tok)
    if not m:
        return [tok]
    return [tok[:m.start()] + v + tok[m.end():] for v in m.group(0)[1:-1].split(",")]


#: a path a document itself marks as still being written is PENDING, not MISSING.  Project
#: memory, "Check the job, not just the file": missing, unfinished and crashed look identical
#: to `ls`, and S29 declared a COMPLETED 126-row result unresolved on a mid-write read.
PENDING_WORDS = ("running", "in flight", "in-flight", "pending", "queued", "launched",
                 "will write", "to be written", "not yet")


def path_audit(docs=DOCS, quiet=False):
    seen, missing, pending = {}, [], []
    for src in docs:
        if not os.path.exists(src):
            continue
        for line in io.open(src, encoding="utf-8").read().splitlines():
            marked = any(w in line.lower() for w in PENDING_WORDS)
            for tok in PATH_RE.findall(line):
                for p in expand(tok):
                    p = p.replace("\\", "/")
                    if p.split("/")[0] in SPRINT_DIRS:
                        s, m = seen.setdefault(p, (set(), [False]))
                        s.add(src)
                        m[0] = m[0] or marked
    for p in sorted(seen):
        srcs, marked = seen[p]
        good = bool(glob.glob(p)) if "*" in p else os.path.exists(p)
        if good:
            OK.append(("path exists: %s" % p, "exists", "yes"))
        elif marked[0]:
            pending.append(p)
            if not quiet:
                print("  pending      %-56s  named as in-flight in %s" % (p, ", ".join(sorted(srcs))))
        else:
            missing.append(p)
            BAD.append(("path claimed but ABSENT: %s" % p, "exists", "*** MISSING ***"))
            if not quiet:
                print("  *** MISSING *** %-56s  named in %s" % (p, ", ".join(sorted(srcs))))
    if not quiet:
        print("  %d distinct paths named across %d S32 documents; %d exist, %d pending, %d MISSING"
              % (len(seen), sum(os.path.exists(d) for d in docs),
                 len(seen) - len(missing) - len(pending), len(pending), len(missing)))
    return seen, missing


if not selftest_only:
    path_audit()

if paths_only:
    print("=" * 98)
    sys.exit(1 if BAD else 0)


# ================================================ AUDIT 1: CROSS-BASIS (carried from S31)
#: every s29 O-ladder item is persisted on BOTH bases, so a claimed value can be matched
#: against both and the DECLARED basis checked against the fit.
_BB = {}
for _r in jsonl("s29/results/s29_O_chain_rows*.jsonl"):
    if _r.get("item") is not None and _r.get("rmsd_chain") == _r.get("rmsd_chain"):
        _BB.setdefault(_r["item"], [[], []])
        _BB[_r["item"]][0].append(float(_r["rmsd_chain"]))
        _BB[_r["item"]][1].append(float(_r.get("rmsd_cloud", float("nan"))))
ITEMS = {k: (mean(v[0]), mean(v[1])) for k, v in _BB.items() if len(v[0]) == 126}
_OTHER = {"chain": "cloud", "cloud": "chain"}


def cross_basis_hit(claimed, basis, tol=5e-4):
    """The ladder item whose OTHER basis `claimed` fits better than the basis it names."""
    o = _OTHER.get(basis)
    if o is None:
        return None
    for it, (ch, cd) in ITEMS.items():
        dec, oth = (ch, cd) if basis == "chain" else (cd, ch)
        if dec != dec or oth != oth:
            continue
        if abs(claimed - oth) <= tol and abs(claimed - oth) < abs(claimed - dec) - 1e-9:
            return it
    return None


# ================================================ AUDIT 2: OBJECT (new, S32-D1)
def object_hit(claimed, obj_name, tol=5e-4):
    """The object `claimed` actually is, when that is NOT the object it is attributed to.

    Catches a number that names the right BASIS ('chain') and the wrong OBJECT -- the S32-L4
    failure mode, where five chain-like means over the same targets span 0.031 A.
    """
    if obj_name is None:
        return None
    dec = OBJECTS.get(obj_name, (None,))[0]
    best, bestd = None, None
    for nm, (v, _s, _d) in OBJECTS.items():
        if v is None or nm == obj_name:
            continue
        d = abs(claimed - v)
        if d <= tol and (bestd is None or d < bestd):
            best, bestd = nm, d
    if best is None:
        return None
    if dec is None:
        return best
    return best if bestd < abs(claimed - dec) - 1e-9 else None


if not selftest_only:
    print()
    print("--- AUDIT 1: does each number MATCH the basis it NAMES? (carried from S31-V D3) ---")
    n1 = 0
    for lbl, cl, bs, _ob in CHECKS:
        it = cross_basis_hit(cl, bs)
        if it:
            FLAG.append(("BASIS MISSTATED: %s is the %s of '%s'" % (cl, _OTHER[bs], it), lbl))
            print("  FLAG %-44s claims %s=%.4f but that is the %s of '%s'"
                  % (lbl[:44], bs, cl, _OTHER[bs], it))
            n1 += 1
    print("  %d ladder items carry both bases; %d cross-basis misstatement(s)" % (len(ITEMS), n1))

    print()
    print("--- AUDIT 2: does each CHAIN number match the OBJECT it is attributed to? (S32-D1) ---")
    n2 = 0
    for lbl, cl, bs, ob in CHECKS:
        hit = object_hit(cl, ob)
        if hit:
            FLAG.append(("OBJECT MISATTRIBUTED: %s is '%s', not '%s'" % (cl, hit, ob), lbl))
            print("  FLAG %-44s %.4f attributed to '%s' but it is '%s'" % (lbl[:44], cl, ob, hit))
            n2 += 1
    print("  %d registered objects; %d misattribution(s)" % (len(OBJECTS), n2))


# ============================== AUDIT 3: cross-job chain claims must carry the bit-identity check
CROSSJOB_RE = re.compile(
    r"(?im)^(?!\s*(?:>|\||#)).*\b(?:built[ -]?chain|chain RMSD|chain mean|on the chain)\b.*$")
BITWORDS = ("bit-identical", "bit identical", "same job", "in the same job", "bit-for-bit",
            "bit-exact", "float64 bits", "paired in one process", "same process")


def crossjob_audit(docs=DOCS, quiet=False):
    """Flag a sentence that quotes a chain DELTA and names neither same-job pairing nor a
    bit-identity check.  S32 rule 20: 'same cloud value' is not enough; 0.03 A of the chain
    is decided by float64 bits."""
    hits = []
    for src in docs:
        if not os.path.exists(src):
            continue
        for ln, line in enumerate(io.open(src, encoding="utf-8").read().splitlines(), 1):
            if not CROSSJOB_RE.match(line):
                continue
            #: only lines that actually quote a signed DELTA below the cross-job floor
            m = re.findall(r"[-+]\s?0\.0[0-2]\d*", line)
            if not m:
                continue
            if any(w in line.lower() for w in BITWORDS):
                continue
            hits.append((src, ln, line.strip()[:110]))
    if not quiet:
        for s, ln, t in hits:
            FLAG.append(("cross-job chain delta < 0.03 A with no bit-identity check", "%s:%d" % (s, ln)))
            print("  FLAG %s:%d  %s" % (s, ln, t))
        print("  %d chain-delta line(s) below the 0.03 A cross-job floor without a bit-identity "
              "or same-job statement" % len(hits))
    return hits


if not selftest_only:
    print()
    print("--- AUDIT 3: chain deltas below the cross-job floor must cite bit-identity (S32 rule 20) ---")
    crossjob_audit()


# ======= AUDIT 4: an aggregate whose FAIL18 and non-FAIL18 strata DISAGREE IN SIGN (S32-D3)
#: Contract rule 12: FAIL18 is an outcome-defined stratum.  S31 measured the effect size rising
#: monotonically with a stratum's circularity and named the gradient "the stratum's definition
#: doing the work".  S32-L3's +0.1872 is the same shape: FAIL18 +1.4879 (0W/18L), the other 108
#: -0.0296 at 0.30x (NOT A RESULT).  An aggregate whose two strata point OPPOSITE WAYS must never
#: be quoted as a property of the instrument without both strata in the same sentence.
STRATA_REG = {}
_sj = load("s32/results/s32_V_orderstat_strata.json")
if _sj:
    for _arm, _d in _sj.items():
        if not isinstance(_d, dict):
            continue
        STRATA_REG[round(float(_d["all_mean"]), 4)] = dict(
            arm=_arm, all=float(_d["all_mean"]), fail18=float(_d["fail18_mean"]),
            rest=float(_d["rest108_mean"]), gate=_d.get("rest108_gate"))
#: what counts as qualifying the quotation.  Either the stratification is named, or the MEDIAN
#: is quoted beside the mean -- which is contract rule 1's own requirement and is what makes the
#: opposite-signed strata visible in the first place.
STRATA_WORDS = ("fail18", "stratum", "strata", "the 108", "other 108", "outcome-defined",
                "circular", "not a result on", "median", "med ", "see d3", "retracted")
CONTEXT = 4          # non-blank lines either side: a claim qualified in its own block is fine


def strata_audit(docs=DOCS, quiet=False):
    """Flag an aggregate whose two strata disagree in sign, quoted without its stratification
    and without its median, anywhere in a +-CONTEXT non-blank-line window."""
    hits = []
    for src in docs:
        if not os.path.exists(src):
            continue
        raw = io.open(src, encoding="utf-8").read().splitlines()
        idx = [k for k, l in enumerate(raw) if l.strip()]          # non-blank line numbers
        pos = {k: i for i, k in enumerate(idx)}
        for ln0, line in enumerate(raw):
            if not line.strip():
                continue
            i = pos[ln0]
            win = " ".join(raw[k] for k in idx[max(0, i - CONTEXT):i + CONTEXT + 1]).lower()
            if any(w in win for w in STRATA_WORDS):
                continue
            ex, why = retraction_exempt(raw, ln0, STRATA_WORDS)
            if ex:
                EXEMPT.append((src, ln0 + 1, "AUDIT 4", why))
                continue
            for val, rec in STRATA_REG.items():
                if rec["fail18"] * rec["rest"] >= 0:      # strata agree: nothing to flag
                    continue
                if ("%.4f" % abs(val)) in line:
                    hits.append((src, ln0 + 1, val, rec, line.strip()[:100]))
                    break
    if not quiet:
        for s, ln, val, rec, t in hits:
            FLAG.append(("aggregate %+.4f has FAIL18 %+.4f vs rest %+.4f -- OPPOSITE SIGNS"
                         % (val, rec["fail18"], rec["rest"]), "%s:%d" % (s, ln)))
            print("  FLAG %s:%d  quotes %+.4f ('%s') without its strata: FAIL18 %+.4f, other 108 "
                  "%+.4f (%s)" % (s, ln, val, rec["arm"], rec["fail18"], rec["rest"], rec["gate"]))
        print("  %d stratified arm(s) registered, %d with opposite-signed strata; %d unqualified "
              "quotation(s)"
              % (len(STRATA_REG), sum(1 for r in STRATA_REG.values() if r["fail18"] * r["rest"] < 0),
                 len(hits)))
    return hits


if not selftest_only:
    print()
    print("--- AUDIT 4: an aggregate whose FAIL18 and non-FAIL18 strata disagree in sign (S32-D3) ---")
    strata_audit()

    print()
    print("--- lane V, S32-D3: the stratification of the order-statistic contrast ---")
    if _sj:
        for _arm, _d in _sj.items():
            if not isinstance(_d, dict):
                continue
            show("  %s" % _arm,
                 "ALL %+0.4f | F18 %+0.4f | 108 %+0.4f" % (_d["all_mean"], _d["fail18_mean"],
                                                           _d["rest108_mean"]),
                 "on the 108: %.2fx -> %s" % (dig(_d, "rest108_compare", "effect_over_mde") or 0,
                                              _d.get("rest108_gate")))
        check("500->128 on the 108 (NOT the outcome-defined stratum)", -0.0296,
              dig(_sj, "500->128", "rest108_mean"), basis="member")
        check("500->128 on FAIL18 (outcome-defined, circular)", 1.4879,
              dig(_sj, "500->128", "fail18_mean"), tol=1e-3, basis="member")
        check("128->75 on the 108", 0.0210, dig(_sj, "128->75", "rest108_mean"), basis="member")
        check("128->75 on FAIL18", 0.3611, dig(_sj, "128->75", "fail18_mean"), tol=1e-3, basis="member")
    gt = load("s32/results/s32_V_orderstat_gate.json")
    if gt:
        k = "SCORE top-128 vs RANDOM 128 of 500"
        check("500->128 MEDIAN effect (opposite sign to the mean)", -0.0380,
              dig(gt, k, "compare", "median_effect"), basis="member")
        check("500->128 n_better (the SCORE wins on this many)", 72,
              dig(gt, k, "compare", "n_better"), tol=0, basis="none")
        show("  rule 10: fraction of SINGLE random draws clearing MDE",
             "%.1f%%" % (100 * (dig(gt, k, "draw_aggregate", "frac_single_draws_clearing_mde") or 0)),
             "the 128->75 arm is 53.6%% -- a coin flip")
        show("  DIS rank of the oracle-best member",
             "mean %.1f / median %.1f" % (dig(gt, "rank_of_oracle_best_member", "mean") or 0,
                                          dig(gt, "rank_of_oracle_best_member", "median") or 0),
             "two moments, 36 apart -- name which")


# ================================================================== SELF-TESTS (rule 5)
print()
print("--- SELF-TESTS: each audit is fed the real historical defect and must TRIP on it ---")
ST = []


def st(name, condition, detail=""):
    ST.append((name, bool(condition)))
    print("  %-62s %s %s" % (name, "PASS" if condition else "*** FAIL -- AUDIT IS BROKEN ***", detail))


# ST1 -- S31-V D3: the STATE headline called the BUILT CHAIN 2.1435 a "CA cloud" number.
_bad = cross_basis_hit(2.1435, "cloud")
_good = cross_basis_hit(2.1458, "cloud")
st("ST1a  chain 2.1435 declared 'cloud' is CAUGHT", _bad is not None, "via '%s'" % _bad)
st("ST1b  cloud 2.1458 declared 'cloud' is CLEAN", _good is None,
   "" if _good is None else "FALSE POSITIVE via '%s'" % _good)

# ST2 -- S32-L4/D1: 3.2105 attributed to a cached scalar.  `rmsd_arm` is the cache's own
# emitted chain (3.2148) and is the object a reader naturally assumes the endpoint is.
_bad2 = object_hit(3.214765, "chain CANONICAL (re-projection)")
_good2 = object_hit(3.210534, "chain CANONICAL (re-projection)")
st("ST2a  the EMITTED chain 3.2148 claimed as the canonical endpoint is CAUGHT",
   _bad2 is not None, "via '%s'" % _bad2)
st("ST2b  the true canonical 3.2105 as the canonical endpoint is CLEAN", _good2 is None,
   "" if _good2 is None else "FALSE POSITIVE via '%s'" % _good2)
_bad3 = object_hit(3.204076, "chain CANONICAL (re-projection)")
st("ST2c  the lam=0 arm 3.2041 claimed as the endpoint is CAUGHT", _bad3 is not None,
   "via '%s'" % _bad3)
_bad4 = object_hit(3.048338, "chain CANONICAL (re-projection)")
st("ST2d  the CLOUD 3.0483 claimed as the endpoint is CAUGHT", _bad4 is not None,
   "via '%s'" % _bad4)

# ST3 -- the path audit must notice an absent path, and not cry wolf about a present one.
_tmp = os.path.join("s32", "results", "_s32_verify_selftest_doc.md")
io.open(_tmp, "w", encoding="utf-8").write(
    "cites `s32/s32_verify.py` and `s32/results/NO_SUCH_FILE_9f2a.json`\n")
try:
    _seen, _miss = path_audit([_tmp], quiet=True)
    st("ST3a  a named-but-absent path is CAUGHT",
       "s32/results/NO_SUCH_FILE_9f2a.json" in _miss)
    st("ST3b  a named-and-present path is CLEAN", "s32/s32_verify.py" not in _miss)
finally:
    # the two BAD rows the self-test just injected are its own, not the sprint's
    BAD[:] = [b for b in BAD if "NO_SUCH_FILE_9f2a" not in str(b)]
    OK[:] = [o for o in OK if "_s32_verify_selftest_doc" not in str(o)]
    os.remove(_tmp)

# ST4 -- the cross-job audit must trip on a real unsupported claim and not on a supported one.
_tmp2 = os.path.join("s32", "results", "_s32_verify_selftest_doc2.md")
io.open(_tmp2, "w", encoding="utf-8").write(
    "The built chain improves by -0.012 A over production.\n"
    "The built chain improves by -0.012 A over production, both projected in the same job.\n")
try:
    _h = crossjob_audit([_tmp2], quiet=True)
    st("ST4a  an unsupported cross-job chain delta of -0.012 is CAUGHT", len(_h) == 1)
    st("ST4b  the same delta WITH a same-job statement is CLEAN", len(_h) == 1,
       "(2 candidate lines, 1 flagged)")
finally:
    os.remove(_tmp2)

# ST6 -- S32-D3: the strata audit must trip on an unqualified quotation of an aggregate whose
# FAIL18 and non-FAIL18 halves point opposite ways, and NOT on the same number stated with its
# stratification.  This is the defect that made S32-L3's headline wrong.
_tmp3 = os.path.join("s32", "results", "_s32_verify_selftest_doc3.md")
io.open(_tmp3, "w", encoding="utf-8").write(
    "The score is worse than random: +0.1872 at 1.06x MDE, WORSE.\n"
    + "filler\n" * (2 * CONTEXT + 2) +
    "The score is +0.1872 overall, but on the 108 non-FAIL18 targets it is -0.0296, NOT A RESULT.\n")
try:
    _h3 = strata_audit([_tmp3], quiet=True)
    st("ST6a  an unqualified +0.1872 is CAUGHT", len(_h3) == 1,
       "" if STRATA_REG else "(no strata registered -- artefact absent)")
    st("ST6b  the same number WITH its strata named is CLEAN", len(_h3) == 1,
       "(2 candidate lines, 1 flagged)")
finally:
    os.remove(_tmp3)

# ST5 -- rule 5 in general: the tie check S31 ran on random floats.  Assert that the data the
# tie claim is made on CAN exhibit a tie.  A tie check on continuous scores never fires.
_s4 = load("s32/results/s32_V_step4_endpoint.json")
_ties = max([r.get("n_tied_at_cut", 0) for r in (dig(_s4, "per_target", default={}) or {}).values()] or [0])
st("ST5   the tie-at-the-cut claim is made on data that CAN tie", _ties > 1,
   "max %d candidates tied at the 75th score" % _ties)

# ST7 -- the RETRACTION EXEMPTION, which is the most dangerous check here because it is the one
# that makes other checks stop firing.  Contract rule 13 retracts in place, so a phrase audit
# would flag every correctly-handled retraction forever until a real flag is invisible.  The
# exemption must therefore satisfy THREE things, not one:
#   (a) a struck claim WITH its replacement stated in the same block is exempt;
#   (b) the SAME claim asserted live below the struck block is still CAUGHT -- otherwise the
#       exemption is a laundering route ("write 'corrected' nearby and say anything");
#   (c) a strike with NO replacement stated is NOT exempt -- the exemption is a positive
#       obligation, and "struck without saying what replaces it" is itself a defect this
#       project has committed before.
_doc = ["> ~~Gain exactly 1 means no noise suppression.~~ STRUCK: gain is 1 inside the active",
        "> hull and exactly 0 outside it -- total suppression, not none.",
        "",
        "Gain exactly 1 means no noise suppression.",
        "",
        "> ~~Gain exactly 1 means no noise suppression.~~ RETRACTED.",
        ""]
_a = retraction_exempt(_doc, 0, A5_CORRECTION)
_b = retraction_exempt(_doc, 3, A5_CORRECTION)
_c = retraction_exempt(_doc, 5, A5_CORRECTION)
st("ST7a  a struck claim WITH its replacement stated is EXEMPT", _a[0] is True, _a[1])
st("ST7b  the same claim asserted LIVE below the block is still CAUGHT", _b[0] is False, _b[1])
st("ST7c  a strike with NO replacement stated is NOT exempt", _c[0] is False, _c[1])


# ======================= AUDIT 8: rows vs DISTINCT pdbs in every results jsonl (standing check)
# Lane P shipped 108 rows over 84 distinct pdbs because surplus launchers duplicated work; rows
# looked like progress.  `len(set(pdbs)) == 126` is the cheap assertion that catches the class.
# A file with a legitimate repeat key (item / draw / shard / variant) is allowed many rows per
# target; a file WITHOUT one and with repeats is flagged.
print()
print("--- AUDIT 8: rows vs distinct pdbs, aggregated over shards (the lane-P duplication class) ---")
REPEAT_KEYS = ("item", "draw", "shard", "variant", "tag", "arm", "k", "s", "eps", "rep", "seed",
               "split", "fam", "start", "scorer", "m", "branch")
#: shards of one logical result are ONE object; the check is on the GROUP, because a shard
#: legitimately holds a slice of the targets and only the union must reach 126.
_SHARD_RE = re.compile(r"(_shard\d+|\.s\d+of\d+|_\d+_\d+)(?=\.jsonl$)")
_groups = {}
for _f in sorted(glob.glob(os.path.join("s32", "results", "*.jsonl"))):
    _g = _SHARD_RE.sub("", _f)
    _groups.setdefault(_g, []).append(_f)
_n8 = 0
for _g, _fs in sorted(_groups.items()):
    _pdbs, _rows, _keys, _perfile = [], 0, set(), []
    for _f in _fs:
        _c = 0
        for _ln in io.open(_f, encoding="utf-8"):
            _ln = _ln.strip()
            if not _ln:
                continue
            try:
                _r = json.loads(_ln)
            except Exception:                                        # noqa: BLE001
                continue
            _rows += 1; _c += 1
            _keys |= set(_r.keys())
            if "pdb" in _r:
                _pdbs.append(_r["pdb"])
        _perfile.append(_c)
    if not _pdbs:
        continue
    _d = len(set(_pdbs))
    _rk = sorted(_keys & set(REPEAT_KEYS))
    #: a repeat key legitimises many rows per target; without one, a repeated pdb is duplicated work
    _dup = 0
    if not _rk:
        _dup = len(_pdbs) - _d
    _bad = _dup > 0
    _short = (_d > 0) and (_rows >= 126) and (_d < 126) and not _rk
    _tag = ""
    if _bad:
        _tag = "   *** %d duplicate pdb row(s), no repeat key ***" % _dup
    elif _short:
        _tag = "   *** %d rows but only %d targets ***" % (_rows, _d)
    print("  %-44s files %d  rows %5d  distinct pdbs %3d  %s%s"
          % (os.path.basename(_g)[:44], len(_fs), _rows, _d,
             ("repeat key %s" % _rk) if _rk else "no repeat key", _tag))
    if _bad or _short:
        FLAG.append(("rows %d over only %d distinct pdbs, repeat keys %s"
                     % (_rows, _d, _rk or "none"), _g))
        _n8 += 1
print("  %d result GROUP(s) with duplicated targets and no repeat key" % _n8)


ST_FAIL = [n for n, ok in ST if not ok]
for n in ST_FAIL:
    FLAG.append(("SELF-TEST FAILED -- the audit it guards is decoration", n))

print()
print("=" * 98)
print("MATCHED: %d    MISMATCHED: %d    KEYS/FILES NOT FOUND: %d    FLAGGED: %d    SELF-TESTS: %d/%d"
      % (len(OK), len(BAD), len(MISSING), len(FLAG), len(ST) - len(ST_FAIL), len(ST)))
for lbl, c, a in BAD:
    print("  MISMATCH  %s: recorded %s vs artefact %s" % (lbl, c, a))
for m in MISSING[:25]:
    print("  NOT FOUND %s" % m)
for why, lbl in FLAG:
    print("  FLAG  %-56s %s" % (why[:56], lbl))
if EXEMPT:
    print("  -- %d struck claim(s) exempted (rule 13, replacement stated in the same block):" % len(EXEMPT))
    for f, ln, aud, why in EXEMPT:
        print("     exempt  %-8s %s:%d  %s" % (aud, f, ln, why))
print("=" * 98)
sys.exit(1 if (BAD or MISSING or ST_FAIL) else 0)
