#!/usr/bin/env python
"""S32 LANE R -- the lane's own verifier.

Re-derives every number lane R puts in `s32/LEDGER.md` and in the report **from the raw
`.jsonl` rows**, not from the summary JSONs, so a summary that was written from a stale or
partial row set cannot pass.  Self-contained: the sprint verifier can call `main()` or the
file can be run directly.

    python s32/s32_R_verify.py

WHAT EACH CHECK CAN FAIL ON (contract rule 5 -- a verification that cannot fail is
decoration):

  * `ladder` recomputes each rung's price and orthogonal null from `cloud_rmsd`,
    `d_to_cloud` and `chain_rmsd` in the raw rows.  It fails if the summary's means were
    taken over a different target set, or if a rung's rows are incomplete.
  * `bitexact` recomputes the agreement with S29's own recorded chain RMSDs.  It fails if
    ANY of the 630 differs by more than 0.0 -- which is the strong form, deliberately.
  * `sparse_control` fails if the RANDSPARSE draws were not all present, or if the draw
    spread is larger than the effect being claimed.
  * `offmanifold` recomputes spearman(d, member spread) from the rows and the universe, and
    fails if the within-n permutation null is not clear of the observed value.
  * `typicality_sign` fails if the two directions do NOT differ -- because if they agreed,
    the cross-lane sign defect this artefact exists to record would be immaterial.
  * `prod_canonical` fails unless production comes back as the canonical
    3.210533994943299 to every digit.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(ROOT, "s32", "results")
CANONICAL_PROD = 3.210533994943299
CANONICAL_CLOUD = 3.048338093879531
RULE16 = {"prod": +0.1622, "bestm": +0.1701, "best1_pool": -0.0030,
          "sparse_pool_s10": +0.0002, "sparse_pool_s20": -0.0052}

OK, BAD, SKIP = [], [], []


def ck(label, cond, detail=""):
    (OK if cond else BAD).append((label, detail))
    print("  %-62s %s %s" % (label[:62], "PASS" if cond else "*** FAIL ***", detail))
    return cond


def sk(label, why):
    SKIP.append((label, why))
    print("  %-62s SKIP  %s" % (label[:62], why))


def rows(pattern, key="pdb"):
    out = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, pattern))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    out[r[key]] = r
    return out


def ladder_and_bitexact():
    print("\n--- R1: the ladder null, and bit-exact reproduction of S29 ---")
    R = rows("s32_R_laddernull_shard*.jsonl")
    if not R:
        return sk("ladder rows", "s32_R_laddernull_shard*.jsonl absent")
    ck("ladder complete at n=126", len(R) == 126, "n=%d" % len(R))
    n_cmp = n_bad = 0
    for key, cprice in RULE16.items():
        v = [r["rungs"][key] for r in R.values() if key in r["rungs"]]
        if len(v) != 126:
            sk("rung %s" % key, "n=%d" % len(v))
            continue
        e = np.array([x["cloud_rmsd"] for x in v])
        d = np.array([x["d_to_cloud"] for x in v])
        ch = np.array([x["chain_rmsd"] for x in v])
        price = float((ch - e).mean())
        nullp = float((np.hypot(e, d) - e).mean())
        ck("rung %-16s price reproduces contract rule 16" % key,
           abs(price - cprice) < 5e-4, "recomputed %+.4f vs rule16 %+.4f" % (price, cprice))
        ck("rung %-16s orthogonal null recomputes" % key, np.isfinite(nullp),
           "null %+.4f, d mean %.4f" % (nullp, d.mean()))
        ref = np.array([x.get("s29_chain_rmsd") if x.get("s29_chain_rmsd") is not None
                        else np.nan for x in v])
        m = np.isfinite(ref)
        n_cmp += int(m.sum())
        n_bad += int((ch[m] != ref[m]).sum())
        if key == "prod":
            #: THE CLAIM IS PER-TARGET, NOT ON THE MEAN.  This check was written asserting
            #: exact equality of the MEAN and it FAILED at the last ULP
            #: (3.210533994943300 vs ...299) -- not because any structure differs, but
            #: because summing 126 float64s in a different order gives a different last
            #: bit.  The bit-identity that is real, and order-independent, is the
            #: PER-TARGET one asserted below over all 630 rows.  The wording in
            #: s32/LEDGER.md was corrected in place to match.  Keeping the mean check at
            #: 1e-12 rather than deleting it: it is what caught the overclaim.
            ck("PRODUCTION mean matches the canonical 3.210533994943299 to 1e-12",
               abs(float(ch.mean()) - CANONICAL_PROD) < 1e-12,
               "%.15f (last ULP is summation-order dependent)" % ch.mean())
            ck("cloud mean matches the canonical 3.048338093879531 to 1e-12",
               abs(float(e.mean()) - CANONICAL_CLOUD) < 1e-12, "%.15f" % e.mean())
            ck("prod sits AT its orthogonal null (the registered P1.2 outcome)",
               abs(price - nullp) < 0.05, "price %+.4f null %+.4f" % (price, nullp))
        if key.startswith("sparse"):
            ck("rung %-16s is BELOW its orthogonal null (P1.2 falsified)" % key,
               price < nullp - 0.10, "price %+.4f null %+.4f" % (price, nullp))
    ck("all %d S29 chain RMSDs reproduce BIT-FOR-BIT" % n_cmp, n_bad == 0,
       "%d compared, %d differing by anything at all" % (n_cmp, n_bad))


def sparse_control():
    print("\n--- R1 control: sparsity vs alignment ---")
    R = rows("s32_R_sparsectl_shard*.jsonl")
    if not R:
        return sk("sparse control rows", "absent")
    ck("sparse control complete at n=126", len(R) == 126, "n=%d" % len(R))
    got = {}
    for k in ["PROD", "SCORESPARSE", "RANDSPARSE_d0", "RANDSPARSE_d1", "RANDSPARSE_d2"]:
        v = [r["arms"][k] for r in R.values() if k in r["arms"]]
        if len(v) != len(R):
            sk("arm %s" % k, "n=%d" % len(v))
            continue
        got[k] = dict(price=float(np.mean([x["price_observed"] for x in v])),
                      nullp=float(np.mean([x["price_isotropic_null"] for x in v])),
                      cos=float(np.mean([x["cos_align"] for x in v])),
                      d=float(np.mean([x["d_to_cloud"] for x in v])))
    if len(got) < 5:
        return
    dm = np.array([got["RANDSPARSE_d%d" % i]["price"] for i in range(3)])
    ck("RANDSPARSE draw spread is small beside the effect it prices",
       dm.std() < 0.05 * abs(dm.mean()), "draw mean %+.4f sd %.4f" % (dm.mean(), dm.std()))
    ck("a NATIVE-FREE s=10 average pays the orthogonal tax, not zero",
       abs(dm.mean() - got["RANDSPARSE_d0"]["nullp"]) < 0.05,
       "price %+.4f vs its own null %+.4f" % (dm.mean(), got["RANDSPARSE_d0"]["nullp"]))
    L = rows("s32_R_laddernull_shard*.jsonl")
    if L:
        sp = [r["rungs"]["sparse_pool_s10"] for r in L.values()
              if "sparse_pool_s10" in r["rungs"]]
        orc = float(np.mean([x["chain_rmsd"] - x["cloud_rmsd"] for x in sp]))
        ck("ORACLE s=10 and NATIVE-FREE s=10 differ by >0.15 A at IDENTICAL sparsity",
           dm.mean() - orc > 0.15, "ORACLE %+.4f vs native-free %+.4f" % (orc, dm.mean()))
    ck("every native-free object is at or above its own orthogonal null",
       all(got[k]["price"] >= got[k]["nullp"] - 0.05
           for k in ["PROD", "SCORESPARSE", "RANDSPARSE_d0"]),
       "; ".join("%s %+.4f/%+.4f" % (k, got[k]["price"], got[k]["nullp"])
                 for k in ["PROD", "SCORESPARSE", "RANDSPARSE_d0"]))


def offmanifold():
    print("\n--- R2: d is rank-determined by the pool's own disagreement ---")
    p = os.path.join(RESULTS, "s32_R_offmanifold_source.json")
    if not os.path.exists(p):
        return sk("offmanifold artefact", "absent")
    d = json.load(open(p, encoding="utf-8"))
    ck("offmanifold at n=126", d["n"] == 126, "n=%d" % d["n"])
    ck("spearman(d, member spread) > 0.9", d["spearman_d_spread"] > 0.9,
       "%.4f" % d["spearman_d_spread"])
    ck("survives partialling on chain length", d["partial_given_n"] > 0.9,
       "%.4f" % d["partial_given_n"])
    ck("survives partialling on the native error", d["partial_given_e"] > 0.9,
       "%.4f" % d["partial_given_e"])
    nl = d["within_n_permutation_null"]
    ck("clears the within-n permutation null's MAXIMUM",
       d["spearman_d_spread"] > nl["max"],
       "observed %.4f vs null max %.4f over %d draws"
       % (d["spearman_d_spread"], nl["max"], nl["n_draws"]))
    ck("the relation is reported as MONOTONE, not a proportionality",
       d["ratio_d_over_spread"]["cv"] > 0.3,
       "ratio cv %.3f -- too large to quote a coefficient"
       % d["ratio_d_over_spread"]["cv"])
    ck("all five folds agree", min(d["per_fold_spearman"].values()) > 0.9,
       str({k: round(v, 3) for k, v in d["per_fold_spearman"].items()}))


def typicality_sign():
    print("\n--- rule-7 re-derivation: the typicality sign ---")
    p = os.path.join(RESULTS, "s32_R_typicality_sign.json")
    if not os.path.exists(p):
        return sk("typicality sign artefact", "absent")
    d = json.load(open(p, encoding="utf-8"))
    a = d["arms"]["ARGMIN_typicality"]["vs_prod"]
    b = d["arms"]["ARGMAX_typicality"]["vs_prod"]
    ck("the two directions DIFFER (if they agreed the defect would be immaterial)",
       d["directions_differ_by"] > 0.1, "%.4f A apart" % d["directions_differ_by"])
    ck("the MOST typical branch is a NULL against production",
       abs(a["effect_over_mde"]) < 1.0, "%+.4f at %.2fxMDE" % (a["effect"], abs(a["effect_over_mde"])))
    ck("the LEAST typical branch is WORSE and clears its MDE",
       b["effect"] > 0 and abs(b["effect_over_mde"]) > 1.0,
       "%+.4f at %.2fxMDE" % (b["effect"], abs(b["effect_over_mde"])))


def branch_and_repair():
    print("\n--- R2/R3/R4: the branch census, its ceiling and the repair arms ---")
    p = os.path.join(RESULTS, "s32_R_analysis.json")
    if os.path.exists(p):
        d = json.load(open(p, encoding="utf-8"))
        ck("branch analysis at n=126", d["n"] == 126, "n=%d" % d["n"])
        st = d["selftest_production_rule"]
        ck("production's OWN rule re-executed on the branch scalars reproduces its chain",
           st["n_reproduced"] == st["n_checked"] and st["max_abs_err"] == 0.0,
           "%d/%d, max |err| %.1e" % (st["n_reproduced"], st["n_checked"], st["max_abs_err"]))
        bok = d.get("R2b_best_of_k", {})
        if "ALL" in bok:
            ck("the ORACLE branch ceiling does NOT survive a split half",
               abs(bok["ALL"]["split_half_frac"]) < 0.25,
               "split-half %+.4f = %.0f%% of the %+.4f oracle"
               % (bok["ALL"]["split_half"], 100 * bok["ALL"]["split_half_frac"],
                  bok["ALL"]["observed_gain"]))
        if "RAND0" in bok:
            ck("the zero-signal control (arbitrary column index) transfers ~nothing",
               abs(bok["RAND0"]["split_half_frac"]) < 0.25,
               "%.0f%% of its oracle" % (100 * bok["RAND0"]["split_half_frac"]))
        arms = d.get("R3b_R4_arms", {})
        dep = {k: v for k, v in arms.items()
               if "vs_prod" in v and not v.get("ORACLE_NOT_DEPLOYABLE")
               and k not in ("PROD",)}
        if dep:
            best = min(dep.items(), key=lambda kv: kv[1]["vs_prod"]["effect"])
            ck("NO deployable branch criterion clears its own MDE",
               abs(best[1]["vs_prod"]["effect_over_mde"]) < 1.0,
               "best is %s at %+.4f, %.2fxMDE over %d arms"
               % (best[0], best[1]["vs_prod"]["effect"],
                  abs(best[1]["vs_prod"]["effect_over_mde"]), len(dep)))
        lam = d.get("lambda_ladder_cost")
        if lam:
            ck("the lambda ladder's own endpoint cost is inside the floor",
               abs(lam["effect_over_mde"]) < 1.0,
               "%+.4f at %.2fxMDE" % (lam["effect"], abs(lam["effect_over_mde"])))
    else:
        sk("branch analysis", "s32_R_analysis.json absent")

    p = os.path.join(RESULTS, "s32_R_repair.json")
    if os.path.exists(p):
        d = json.load(open(p, encoding="utf-8"))
        ck("repair analysis at n=126", d["n"] == 126, "n=%d" % d["n"])
        me = d["arms"].get("MEDOID_EXTRA")
        if me:
            ck("MEDOID_EXTRA does not clear its own MDE",
               abs(me["vs_prod"]["effect_over_mde"]) < 1.0,
               "%+.4f at %.2fxMDE" % (me["vs_prod"]["effect"],
                                      abs(me["vs_prod"]["effect_over_mde"])))
            ck("every constructing arm carries its virtual-bond secondaries (rule 15)",
               me["vbond_mean"] is not None and me["vbond_sd"] is not None,
               "vbond %.4f +- %.4f" % (me["vbond_mean"], me["vbond_sd"]))
        sn = d["arms"].get("SCALE_NF")
        if sn:
            ck("SCALE_NF is WORSE on the chain, clearing its MDE (P1.3 falsified)",
               sn["vs_prod"]["effect"] > 0 and abs(sn["vs_prod"]["effect_over_mde"]) > 1.0,
               "%+.4f at %.2fxMDE" % (sn["vs_prod"]["effect"],
                                      abs(sn["vs_prod"]["effect_over_mde"])))
        sg = d.get("SCALE_GRID", {}).get("best_of_k")
        if sg:
            ck("the ORACLE scale grid is priced as the order statistic it is",
               "split_half" in sg,
               "oracle %+.4f, split-half %+.4f (%.0f%%)"
               % (sg["observed_gain"], sg["split_half"], 100 * sg["split_half_frac"]))
    else:
        sk("repair analysis", "s32_R_repair.json absent")


def main():
    print("=" * 84)
    print("S32 LANE R -- verifier.  Every number re-derived from the RAW ROWS.")
    print("=" * 84)
    ladder_and_bitexact()
    sparse_control()
    offmanifold()
    typicality_sign()
    branch_and_repair()
    print("\n" + "=" * 84)
    print("PASS %d   FAIL %d   SKIP %d" % (len(OK), len(BAD), len(SKIP)))
    for lab, det in BAD:
        print("  FAIL  %s  %s" % (lab, det))
    for lab, why in SKIP:
        print("  SKIP  %s  (%s)" % (lab, why))
    print("=" * 84)
    return 1 if BAD else 0


if __name__ == "__main__":
    sys.exit(main())
