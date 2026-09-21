#!/usr/bin/env python
"""S32 LANE R -- analysis of the branch census.  R2a / R2b / R3a / R3b / R4-from-rows.

PREREG: s32/PREREG_S32_R.md at 02754f5a.  BASIS: built-chain CA-RMSD, tuning126, n=126.
Every arm in here is a SELECTION over branches that were all built in the same job from the
same cloud, so every contrast is paired inside one process.

SELF-TEST (contract rule 5): production's own selection rule is re-executed on the recorded
branch scalars and the resulting chain RMSD is checked against `prod_chain`, which was
produced by `lam_path` itself.  A wrong branch-set definition makes this fail; it is not
decoration.
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

from s24 import stats_lib as ST                                          # noqa: E402
from s12 import instrument as I                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")

#: criteria, all "lower is better" as registered in the PREREG.
#: NOTE ON `typicality`: it is a DISTANCE to the pool (mean CA-RMSD to the top-75), so lower
#: = MORE typical.  Minimising it selects the most typical branch.  Lane V's adversary script
#: maximised it and so measured the LEAST typical branch; both directions are reported in
#: `s32/results/s32_R_typicality_sign.json`.
CRITERIA = ["rama_nlp", "rama20_nlp", "ramah", "posphi_frac", "disto_risk", "disto_mae",
            "legacy", "typicality", "obj1", "obj0", "d_to_C"]
SUBSET_SIZES = (1, 2, 4, 8, 16, 32, 64, 128)

#: THE ANALYSED BRANCH SET, fixed for every target so that no quantity is compared across a
#: changed object.  Rows written before 09:20 on 2026-09-21 also carry RAND1/RAND2; they are
#: dropped here.  158 branches per target: 4 + 4 + 75 + 75.
FAMILIES = ("GEN4", "GEN4D", "MEM75", "RAND0")


STRUCTS = os.path.join(RESULTS, "s32_R_structs")


def components(pdb, thresh=1e-3):
    """CONNECTED COMPONENTS of the `pairwise RMSD < thresh` adjacency, by union-find.

    Recomputed here from the persisted torsions rather than trusted from the rows: rows
    written before 2026-09-21 carry a GREEDY-SEED count, which overwrites labels an earlier
    seed assigned and is therefore not a component count.  Lane V's adversary pass found it.
    Returns (labels, n_components) or (None, None) if the npz is absent.

    SELF-TEST (contract rule 5): this CAN differ from the stored labels, and the analysis
    prints how often it does.  On data where every branch is isolated the two agree, which
    is exactly the input that would have hidden the bug.
    """
    f = os.path.join(STRUCTS, f"{pdb}.npz")
    if not os.path.exists(f):
        return None, None, None
    from core import project as pj
    with np.load(f, allow_pickle=True) as z:
        PH = np.asarray(z["phi"], float)
        PS = np.asarray(z["psi"], float)
        fam = np.asarray(z["fam"]).astype(str)
    keep = np.isin(fam, FAMILIES)
    PH, PS = PH[keep], PS[keep]
    CA = np.asarray(pj.build_ca_exact(PH, PS), float)
    B = len(CA)
    par = np.arange(B)

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for a in range(B):
        row = I.kabsch_rmsd_batch(CA, CA[a])
        for b in np.where(row < thresh)[0]:
            ra, rb = find(a), find(int(b))
            if ra != rb:
                par[ra] = rb
    roots = np.array([find(a) for a in range(B)])
    _, lab = np.unique(roots, return_inverse=True)
    return lab, int(lab.max()) + 1, keep


def load():
    rows = {}
    for p in sorted(glob.glob(os.path.join(RESULTS, "s32_R_branches_shard*.jsonl"))):
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    rows[r["pdb"]] = r
    #: repair the greedy-seed labels in place
    n_changed, n_ok, per_branch = 0, 0, None
    for p, r in rows.items():
        lab, nc, keep = components(p)
        if lab is None:
            continue
        if per_branch is None:
            per_branch = [k for k, v in r.items()
                          if isinstance(v, list) and len(v) == r["n_branches"]]
        for k in per_branch:
            if isinstance(r.get(k), list) and len(r[k]) == r["n_branches"]:
                r[k] = [x for x, m in zip(r[k], keep) if m]
        if nc != r.get("n_distinct"):
            n_changed += 1
        else:
            n_ok += 1
        r["cluster"] = lab.tolist()
        r["n_distinct_greedy_BUGGY"] = r.get("n_distinct")
        r["n_distinct"] = nc
        r["n_branches_all_families"] = r["n_branches"]
        r["n_branches"] = int(keep.sum())
    return rows, {"n_relabelled": n_changed, "n_agreed": n_ok,
                  "note": "union-find components recomputed from persisted torsions; the "
                          "stored count was a greedy-seed count (lane V, 2026-09-21)"}


def argmin_tied_mean(score, outcome, atol=1e-12, rtol=1e-9):
    """Mean outcome over the TIED argmin set.  np.argmin on a tied criterion reads the
    array order, and in S12 that order was the ORACLE sort -- it invented a winner."""
    s = np.asarray(score, float)
    o = np.asarray(outcome, float)
    ok = np.isfinite(s)
    if not ok.any():
        return float("nan"), 0
    s2 = np.where(ok, s, np.inf)
    m = s2.min()
    tie = np.isclose(s2, m, atol=atol, rtol=rtol)
    return float(o[tie].mean()), int(tie.sum())


def main():
    rows, clusterfix = load()
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in rows]
    n = len(pdbs)
    folds = np.array([rows[p]["fold"] for p in pdbs], int)
    out = {"n": n, "pdbs": pdbs, "basis": "built chain CA-RMSD, tuning126",
           "prereg_commit": "02754f5a", "cluster_relabel": clusterfix}
    if n < 126:
        out["INCOMPLETE"] = True

    prod = np.array([rows[p]["prod_chain"] for p in pdbs])
    cloud = np.array([rows[p]["cloud_rmsd"] for p in pdbs])
    prodfit = np.array([rows[p]["prod_fit_chain"] for p in pdbs])
    prodcache = np.array([rows[p]["prod_cachecloud_chain"] for p in pdbs])
    dprodcloud = np.array([rows[p]["prod_d_to_cloud"] for p in pdbs])

    out["headline"] = {
        "PROD_chain_mean": float(prod.mean()), "PROD_chain_se": float(prod.std(ddof=1) / math.sqrt(n)),
        "PROD_chain_median": float(np.median(prod)),
        "PROD_fit_lam0_chain_mean": float(prodfit.mean()),
        "cloud_mean": float(cloud.mean()),
        "price_mean": float((prod - cloud).mean()),
        "price_se": float((prod - cloud).std(ddof=1) / math.sqrt(n)),
        "price_median": float(np.median(prod - cloud)),
        "n_price_positive": int((prod > cloud).sum()),
        "PROD_cachecloud_chain_mean": float(prodcache.mean()),
        "cloud_cache_delta_max": float(max(rows[p]["cloud_cache_delta"] for p in pdbs)),
        "cloud_recon_err_max": float(max(rows[p]["cloud_recon_err"] for p in pdbs)),
    }
    #: the lambda ladder's OWN endpoint cost, both arms from the SAME cloud in the SAME job.
    #: the coordinator quotes +0.0107 from the production cache; 0.0134 is the MEAN chain
    #: floor, so this needs its MDE before anything is read into it.
    out["lambda_ladder_cost"] = ST.compare(prod, prodfit, folds, names=pdbs,
                                           label="lam=0.3 arm - lam=0 arm (same cloud, same job)")

    #: the 1e-14 cloud perturbation, measured as a chain effect
    out["cloud_perturbation_floor"] = ST.compare(prodcache, prod, folds, names=pdbs,
                                                 label="PROD(cache cloud) - PROD(s29 cloud)")
    dd = np.abs(prodcache - prod)
    out["cloud_perturbation_floor"].update(
        per_target_abs_mean=float(dd.mean()), per_target_abs_p90=float(np.percentile(dd, 90)),
        per_target_abs_max=float(dd.max()), n_moved_1e3=int((dd > 1e-3).sum()),
        n_moved_0p1=int((dd > 0.1).sum()))

    # ---------------------------------------------------------------- R1 isotropic null
    #: if the projection displacement were orthogonal to the cloud's native error the chain
    #: error would be sqrt(e^2 + d^2).  Matched to the operator's OWN displacement.
    null = np.sqrt(cloud ** 2 + dprodcloud ** 2)
    out["R1_isotropic_null"] = {
        "d_to_cloud_mean": float(dprodcloud.mean()), "d_to_cloud_median": float(np.median(dprodcloud)),
        "null_chain_mean": float(null.mean()),
        "null_price_mean": float((null - cloud).mean()),
        "observed_price_mean": float((prod - cloud).mean()),
        "compare_PROD_vs_null": ST.compare(prod, null, folds, names=pdbs,
                                           label="PROD chain - isotropic null"),
    }

    # ------------------------------------------------- R2 the ball bound on the ceiling
    #: Kabsch RMSD is a metric on shape space, so every branch X satisfies
    #:     |RMSD(X, nat) - RMSD(C, nat)| <= RMSD(X, C).
    #: Every branch is a near-optimal projection, so RMSD(X, C) ~ d.  The BEST a branch can
    #: possibly be is therefore cloud_rmsd - d, and that requires the projection displacement
    #: to point exactly at the native -- a measure-zero direction.  This is the hard bound on
    #: everything R2 can win, and it is free.
    dbr = np.array([float(np.median(rows[p]["d_to_C"])) for p in pdbs])
    ball = np.maximum(cloud - dbr, 0.0)
    out["R2_ball_bound"] = {
        "branch_d_to_C_median_of_medians": float(np.median(dbr)),
        "branch_d_to_C_mean": float(dbr.mean()),
        "ball_bound_mean": float(ball.mean()),
        "ball_bound_note": "cloud_rmsd - d; a LOWER BOUND on any branch, not achievable",
        "isotropic_null_mean": float(np.sqrt(cloud ** 2 + dbr ** 2).mean()),
    }

    # ---------------------------------------------------------------- R2a the census
    cen = []
    for p in pdbs:
        r = rows[p]
        rm = np.array(r["rmsd_nat"]); cl = np.array(r["cluster"]); ob = np.array(r["obj1"])
        reps = np.array([np.where(cl == c)[0][0] for c in np.unique(cl)])
        rr, oo = rm[reps], ob[reps]
        o = np.argsort(oo)
        cen.append(dict(pdb=p, n=r["n"], n_branches=r["n_branches"], n_distinct=r["n_distinct"],
                        rmsd_min=float(rr.min()), rmsd_max=float(rr.max()),
                        rmsd_spread=float(rr.max() - rr.min()),
                        rmsd_iqr=float(np.percentile(rr, 75) - np.percentile(rr, 25)),
                        rmsd_p90_p10=float(np.percentile(rr, 90) - np.percentile(rr, 10)),
                        obj_gap_best_runnerup=float(oo[o[1]] - oo[o[0]]) if len(o) > 1 else float("nan"),
                        obj_spread=float(oo.max() - oo.min()),
                        prod_chain=float(r["prod_chain"])))
    sp = np.array([c["rmsd_spread"] for c in cen])
    nd = np.array([c["n_distinct"] for c in cen])
    gp = np.array([c["obj_gap_best_runnerup"] for c in cen])
    out["R2a_census"] = {
        "rows": cen,
        "n_distinct_median": float(np.median(nd)), "n_distinct_min": int(nd.min()),
        "n_targets_ge2_distinct": int((nd >= 2).sum()),
        "rmsd_spread_mean": float(sp.mean()), "rmsd_spread_median": float(np.median(sp)),
        "rmsd_spread_p90": float(np.percentile(sp, 90)), "rmsd_spread_max": float(sp.max()),
        "n_spread_gt_0p05": int((sp > 0.05).sum()), "n_spread_gt_0p3": int((sp > 0.3).sum()),
        "n_spread_gt_1p0": int((sp > 1.0).sum()),
        "obj_gap_median": float(np.nanmedian(gp)),
        "n_obj_gap_under_1e6": int(np.nansum(gp < 1e-6)),
        "n_obj_gap_under_1e3": int(np.nansum(gp < 1e-3)),
    }

    # ---------------------------------------------------------------- self-test rule 5
    st = {"n_checked": 0, "n_reproduced": 0, "max_abs_err": 0.0, "errs": []}
    for p in pdbs:
        r = rows[p]
        fam = np.array(r["fam"]); rm = np.array(r["rmsd_nat"])
        o0 = np.array(r["obj0"]); o1 = np.array(r["obj1"])
        g = np.where(fam == "GEN4")[0]; gd = np.where(fam == "GEN4D")[0]
        if not len(g) or not len(gd):
            continue
        k0 = g[int(np.argmin(o0[g]))]
        cand = np.concatenate([[k0], gd])
        pick = cand[int(np.argmin(o1[cand]))]
        e = abs(float(rm[pick]) - float(r["prod_chain"]))
        st["n_checked"] += 1
        st["n_reproduced"] += int(e < 1e-6)
        st["max_abs_err"] = max(st["max_abs_err"], e)
        if e >= 1e-6:
            st["errs"].append([p, float(rm[pick]), float(r["prod_chain"]), e])
    st["errs"] = st["errs"][:15]
    out["selftest_production_rule"] = st

    # ---------------------------------------------------------------- R2b ORACLE ceiling
    fams = list(FAMILIES)
    orc = {}
    for f in fams + ["ALL", "GEN4+GEN4D", "GEN4+MEM75"]:
        v = []
        for p in pdbs:
            r = rows[p]
            fam = np.array(r["fam"]); rm = np.array(r["rmsd_nat"])
            if f == "ALL":
                m = np.ones(len(fam), bool)
            elif "+" in f:
                m = np.isin(fam, f.split("+"))
            else:
                m = fam == f
            v.append(rm[m].min())
        v = np.array(v)
        orc[f] = {"ORACLE_NOT_DEPLOYABLE": True, "mean": float(v.mean()),
                  "se": float(v.std(ddof=1) / math.sqrt(n)), "median": float(np.median(v)),
                  "vs_prod": ST.compare(v, prod, folds, names=pdbs,
                                        label="ORACLE best branch (%s) - PROD" % f)}
    out["R2b_oracle_ceiling"] = orc

    #: the order-statistic curve: ORACLE min over a RANDOM SUBSET of K branches.  If this is
    #: still falling at the largest K the ceiling is a draw-size artefact (contract rule 9).
    rng = np.random.default_rng(3222)
    curve = {}
    for K in SUBSET_SIZES:
        draws = []
        for _ in range(24):
            v = []
            for p in pdbs:
                rm = np.array(rows[p]["rmsd_nat"])
                k = min(K, len(rm))
                v.append(rm[rng.choice(len(rm), size=k, replace=False)].min())
            draws.append(float(np.mean(v)))
        curve[str(K)] = {"ORACLE_NOT_DEPLOYABLE": True, "draw_mean": float(np.mean(draws)),
                         "draw_sd": float(np.std(draws)), "n_draws": len(draws)}
    curve["ALL"] = {"ORACLE_NOT_DEPLOYABLE": True,
                    "draw_mean": orc["ALL"]["mean"], "draw_sd": 0.0, "n_draws": 1}
    out["R2b_order_statistic_curve"] = curve

    #: BEST-OF-K PRICING AND SPLIT-HALF TRANSFER (contract rule 9).  A column of `M` is one
    #: SETTING held fixed across targets.  For GEN4/GEN4D a column is one named generic start;
    #: for MEM75 it is "the r-th distogram-ranked member's torsions" -- both are transferable
    #: labels.  For RAND the column index is arbitrary BY CONSTRUCTION, so its split-half
    #: transfer is the built-in zero-signal control and must come out at ~0.
    bok = {}
    for f in list(FAMILIES) + ["ALL"]:
        cols = []
        ok = True
        for p in pdbs:
            r = rows[p]
            fam = np.array(r["fam"]); rm = np.array(r["rmsd_nat"])
            m = np.ones(len(fam), bool) if f == "ALL" else (fam == f)
            cols.append(rm[m])
        w = min(len(c) for c in cols)
        if w < 2:
            continue
        M = np.array([c[:w] for c in cols])
        bok[f] = ST.best_of_k_within(M, n_boot=300, seed_parts=("s32R", f))
        bok[f]["K"] = int(w)
    #: production's own 8-candidate set, the S31-D comparison point
    cols = []
    for p in pdbs:
        r = rows[p]
        fam = np.array(r["fam"]); rm = np.array(r["rmsd_nat"])
        cols.append(np.concatenate([rm[fam == "GEN4"], rm[fam == "GEN4D"]]))
    w = min(len(c) for c in cols)
    bok["PROD8"] = ST.best_of_k_within(np.array([c[:w] for c in cols]), n_boot=300,
                                       seed_parts=("s32R", "PROD8"))
    bok["PROD8"]["K"] = int(w)
    bok["_note"] = ("ORACLE / NOT DEPLOYABLE.  `split_half` is the number to quote: it nulls "
                    "itself.  RAND0's column index is arbitrary by construction and is the "
                    "zero-signal control for the transfer statistic.")
    out["R2b_best_of_k"] = bok

    # ---------------------------------------------------------------- R3a in-band skill
    from scipy.stats import spearmanr
    inband = {}
    for c in CRITERIA:
        rho, nties = [], []
        for p in pdbs:
            r = rows[p]
            cl = np.array(r["cluster"])
            reps = np.array([np.where(cl == q)[0][0] for q in np.unique(cl)])
            x = np.array(r[c], float)[reps]
            y = np.array(r["rmsd_nat"], float)[reps]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 3 or np.ptp(x[ok]) == 0:
                rho.append(np.nan); continue
            rho.append(float(spearmanr(x[ok], y[ok]).statistic))
        rho = np.array(rho, float)
        ok = np.isfinite(rho)
        rr = rho[ok]; ff = folds[ok]
        F = np.array(sorted(set(ff.tolist())))
        rg2 = np.random.default_rng(abs(hash(c)) % 2 ** 31)
        bf = np.array([np.concatenate([rr[ff == q] for q in rg2.choice(F, len(F), replace=True)]).mean()
                       for _ in range(4000)])
        inband[c] = {"mean_rho": float(rr.mean()), "median_rho": float(np.median(rr)),
                     "n": int(ok.sum()), "se": float(rr.std(ddof=1) / math.sqrt(len(rr))),
                     "ci95_fold": [float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))],
                     "n_positive": int((rr > 0).sum()),
                     "per_fold": {int(q): float(rr[ff == q].mean()) for q in F}}
    out["R3a_inband_rho"] = inband

    # ---------------------------------------------------------------- R3b/R4 directional
    arms = {}

    def add(name, vals, oracle=False, note=""):
        v = np.asarray(vals, float)
        arms[name] = {"mean": float(v.mean()), "se": float(v.std(ddof=1) / math.sqrt(n)),
                      "median": float(np.median(v)), "ORACLE_NOT_DEPLOYABLE": bool(oracle),
                      "note": note,
                      "vs_prod": ST.compare(v, prod, folds, names=pdbs,
                                            label="%s - PROD" % name)}
        return v

    add("PROD", prod, note="baseline, lam_path(extra=None), same job")
    for c in CRITERIA:
        for famset, tag in ((None, "ALL"), (["GEN4", "MEM75"], "GEN4+MEM75")):
            v, nt = [], []
            for p in pdbs:
                r = rows[p]
                fam = np.array(r["fam"])
                m = np.ones(len(fam), bool) if famset is None else np.isin(fam, famset)
                a, k = argmin_tied_mean(np.array(r[c], float)[m], np.array(r["rmsd_nat"])[m])
                v.append(a); nt.append(k)
            add("SEL_%s_%s" % (c, tag), v,
                note="argmin over tied set averaged; mean tie size %.2f" % np.mean(nt))

    #: controls
    for p_i, tag in ((None, "ALL"), (["GEN4", "MEM75"], "GEN4+MEM75")):
        v = []
        for p in pdbs:
            r = rows[p]
            fam = np.array(r["fam"])
            m = np.ones(len(fam), bool) if p_i is None else np.isin(fam, p_i)
            v.append(np.array(r["rmsd_nat"])[m].mean())
        add("BRANCHMEAN_%s" % tag, v, note="zero-skill control: the branch set's own mean")

    rngc = np.random.default_rng(99117)
    draws = []
    for d in range(5):
        v = [np.array(rows[p]["rmsd_nat"])[rngc.integers(0, rows[p]["n_branches"])] for p in pdbs]
        draws.append(float(np.mean(v)))
        if d == 0:
            add("RANDBRANCH_draw0", v, note="one random branch, draw 0 of 5")
    arms["RANDBRANCH_distribution"] = {"draw_means": draws, "draw_mean": float(np.mean(draws)),
                                       "draw_sd": float(np.std(draws)), "n_draws": len(draws)}

    #: R4: the medoid-torsion start alone, and the objective rule over wider start sets
    v = []
    for p in pdbs:
        r = rows[p]
        fam = np.array(r["fam"]); tag = np.array(r["tag"])
        k = np.where((fam == "MEM75") & (tag == "mem_%d" % r["medoid_pos"]))[0]
        v.append(float(np.array(r["rmsd_nat"])[k[0]]) if len(k) else np.nan)
    add("MEDOID_ONLY", v, note="the branch reached from the medoid member's own torsions")

    for famset, tag in ((["GEN4", "MEM75"], "GEN4+MEM75"), (None, "ALL"),
                        (["RAND0"], "RAND0")):
        v = []
        for p in pdbs:
            r = rows[p]
            fam = np.array(r["fam"])
            m = np.ones(len(fam), bool) if famset is None else np.isin(fam, famset)
            a, _ = argmin_tied_mean(np.array(r["obj1"], float)[m], np.array(r["rmsd_nat"])[m])
            v.append(a)
        add("OBJARGMIN_%s" % tag, v, note="production's own selection rule over a wider start set")

    #: ORACLE arms, for the ceiling only
    add("ORACLE_BEST_ALL", [np.array(rows[p]["rmsd_nat"]).min() for p in pdbs], oracle=True)
    add("ORACLE_BEST_GEN4", [np.array(rows[p]["rmsd_nat"])[np.array(rows[p]["fam"]) == "GEN4"].min()
                             for p in pdbs], oracle=True)
    out["R3b_R4_arms"] = arms

    ST.save_atomic(os.path.join(RESULTS, "s32_R_analysis.json"), out, module_file=__file__)
    print(json.dumps({k: out[k] for k in ("n", "headline", "R1_isotropic_null",
                                          "selftest_production_rule")},
                     indent=1, default=str)[:4000])
    print("\n== lambda ladder cost (lam=0.3 arm - lam=0 arm, same cloud, same job) ==")
    print(ST.fmt(out["lambda_ladder_cost"]))
    print("\n== 1e-14 cloud perturbation, as a CHAIN effect ==")
    print(json.dumps({k: out["cloud_perturbation_floor"][k] for k in
                      ("per_target_abs_mean", "per_target_abs_p90", "per_target_abs_max",
                       "n_moved_1e3", "n_moved_0p1", "effect", "effect_over_mde")}, indent=1))
    print("\n== R2a census ==")
    print(json.dumps({k: v for k, v in out["R2a_census"].items() if k != "rows"}, indent=1))
    print("\n== R2b order-statistic curve (ORACLE / NOT DEPLOYABLE) ==")
    for k, v in out["R2b_order_statistic_curve"].items():
        print("  K=%-5s %.4f +- %.4f" % (k, v["draw_mean"], v["draw_sd"]))
    print("\n== R2b best-of-K pricing / split-half transfer (ORACLE) ==")
    for k, v in out["R2b_best_of_k"].items():
        if k.startswith("_"):
            continue
        print("  %-7s K=%-4d oracle %+.4f  across-target null %+.4f (%.0f%%)  "
              "SPLIT-HALF %+.4f (%.0f%% of oracle)  k_eff %.1f"
              % (k, v["K"], v["observed_gain"], v["null_across_targets"],
                 100 * v["share_accounted"], v["split_half"], 100 * v["split_half_frac"],
                 v["k_eff"]))
    print("\n== R3a in-band rho (positive = criterion ranks RMSD correctly) ==")
    for c, v in sorted(inband.items(), key=lambda kv: -kv[1]["mean_rho"]):
        print("  %-12s rho %+.4f  median %+.4f  foldCI [%+.4f, %+.4f]  %d/%d positive"
              % (c, v["mean_rho"], v["median_rho"], v["ci95_fold"][0], v["ci95_fold"][1],
                 v["n_positive"], v["n"]))
    print("\n== R3b/R4 arms vs PROD (LOWER IS BETTER) ==")
    for k, v in sorted(arms.items(), key=lambda kv: kv[1].get("mean", 9e9)):
        if "vs_prod" not in v:
            continue
        c = v["vs_prod"]
        print("  %-28s %.4f  eff %+.4f  %.2fxMDE  %dW/%dL  %s%s"
              % (k, v["mean"], c["effect"], abs(c["effect_over_mde"]),
                 c["n_better"], c["n_worse"], c["verdict"][:46],
                 "   [ORACLE]" if v["ORACLE_NOT_DEPLOYABLE"] else ""))
    print("\nwrote", os.path.join(RESULTS, "s32_R_analysis.json"))


if __name__ == "__main__":
    main()
