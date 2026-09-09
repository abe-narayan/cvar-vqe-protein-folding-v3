"""s23/agentA_report.py -- analysis/report over s23/results/agentA.json, per PREREG_A.md.

Produces the full stats package for all three experiments: paired bootstrap iid + fold-clustered
CI, SE, MDE = 2.8016*SE (effect reported as a multiple of its own MDE), W/L, worst-target delta,
mean virtual bond length. Nested (leave-one-fold-out) CV for Experiment 3's fitted parameters,
following `s22/routercv.py`'s convention. Writes `s23/results/agentA_report.json` (all numbers,
machine-readable) and prints the same as text.
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
RES = os.path.join(HERE, "results")

from s15 import seed as SD   # noqa: E402

B = 4000


def _stat(d, folds, rng):
    d = np.asarray(d, float); n = len(d)
    se = float(d.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    mde = 2.8016 * se
    b = d[rng.integers(0, n, size=(B, n))].mean(1)
    ci_iid = (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))
    F = sorted(set(int(f) for f in folds))
    fb = []
    for _ in range(B):
        chosen = rng.choice(F, size=len(F), replace=True)
        fb.append(np.concatenate([d[folds == q] for q in chosen]).mean())
    fb = np.array(fb)
    ci_fold = (float(np.percentile(fb, 2.5)), float(np.percentile(fb, 97.5)))
    w = int((d < 0).sum()); l = int((d > 0).sum())
    return dict(n=n, mean=float(d.mean()), se=se, mde=mde, ci_iid=ci_iid, ci_fold=ci_fold,
                w=w, l=l, worst=float(d.max()), best=float(d.min()),
                mde_multiple=float(abs(d.mean()) / mde) if mde > 0 else float("nan"))


def _fmt(s):
    return ("%+.4f SE %.4f MDE %.4f (%.2fx)  iid[%+.3f,%+.3f] fold[%+.3f,%+.3f]  %dW/%dL  worst %+.3f"
            % (s["mean"], s["se"], s["mde"], s["mde_multiple"], s["ci_iid"][0], s["ci_iid"][1],
               s["ci_fold"][0], s["ci_fold"][1], s["w"], s["l"], s["worst"]))


def load():
    d = json.load(open(os.path.join(RES, "agentA.json")))
    assert d["complete"], "agentA.json is not marked complete"
    return d


def report():
    d = load()
    rows = d["rows"]
    folds = np.array([r["fold"] for r in rows], int)
    inc = np.array([r["avg_75"] for r in rows], float)
    vb_inc = np.array([r["vb_avg_75"] for r in rows], float)
    rng = SD.stable_rng("s23agentA", "report")
    out = {"n": len(rows), "incumbent_mean": float(inc.mean()), "incumbent_vbond": float(vb_inc.mean())}

    print("=" * 100)
    print("WORKSTREAM A -- ENSEMBLE AGGREGATION -- n=%d, incumbent avg_75 = %.4f (vbond %.3f, physical 3.805)"
          % (len(rows), inc.mean(), vb_inc.mean()))
    print("=" * 100)

    # ============================================================= EXPERIMENT 1
    print("\n--- EXPERIMENT 1: CLUSTER-RESTRICTED AVERAGING ---\n")
    exp1_out = {}
    cells = sorted(rows[0]["exp1"].keys())
    for cell in cells:
        print("  [%s]  n_clusters(median)=%s  sizes(median)=%s" % (
            cell, int(np.median([r["exp1"][cell]["n_clusters"] for r in rows])),
            rows[len(rows) // 2]["exp1"][cell]["sizes"]))
        exp1_out[cell] = {}
        rules = ["size", "score", "combined", "random", "oracle"]
        if "random_same_size_partition" in rows[0]["exp1"][cell]:
            rules.append("random_same_size_partition")
        for rule in rules:
            x = np.array([r["exp1"][cell][rule] for r in rows], float)
            s = _stat(x - inc, folds, rng)
            vb = np.mean([r["exp1"][cell].get("vb_" + rule, np.nan) for r in rows]) if rule != "random_same_size_partition" else float("nan")
            flag = ""
            if rule not in ("oracle",):
                if s["ci_iid"][1] < 0 and s["mde_multiple"] > 1.0:
                    flag = "  <<< BEATS INCUMBENT"
            tag = "ORACLE-CEILING" if rule == "oracle" else rule
            print("    %-26s mean=%.4f vbond=%.3f  %s%s" % (tag, x.mean(), vb, _fmt(s), flag))
            exp1_out[cell][rule] = dict(arm_mean=float(x.mean()), vbond=float(vb), **s)
        print()
    out["exp1"] = exp1_out

    # ============================================================= EXPERIMENT 2
    print("\n--- EXPERIMENT 2: ENSEMBLE-WIDTH SWEEP (plain vs cluster-restricted, k=3 combined) ---\n")
    exp2_out = {}
    fail_thresh = d.get("fail_thresh", 5.0)
    print("  %-6s%10s%10s%10s%10s%9s  |  %10s%10s%10s%10s%9s   vs avg_m" % (
        "m", "avg_mean", "avg_med", "avg_worst", "avg_var", "avg_fail",
        "cl_mean", "cl_med", "cl_worst", "cl_var", "cl_fail"))
    for m in d["m_ladder"]:
        key = "m%d" % m
        a = np.array([r["exp2"][key]["avg_m"] for r in rows], float)
        c = np.array([r["exp2"][key]["clust_m"] for r in rows], float)
        s = _stat(c - a, folds, rng)
        row_out = dict(
            avg_mean=float(a.mean()), avg_median=float(np.median(a)), avg_worst=float(a.max()),
            avg_var=float(a.var()), avg_fail_rate=float((a > fail_thresh).mean()),
            clust_mean=float(c.mean()), clust_median=float(np.median(c)), clust_worst=float(c.max()),
            clust_var=float(c.var()), clust_fail_rate=float((c > fail_thresh).mean()),
            clust_vs_avg=s)
        print("  %-6d%10.4f%10.4f%10.4f%10.4f%9.3f  |  %10.4f%10.4f%10.4f%10.4f%9.3f   %+.4f (%.2fx MDE)"
              % (m, a.mean(), np.median(a), a.max(), a.var(), (a > fail_thresh).mean(),
                 c.mean(), np.median(c), c.max(), c.var(), (c > fail_thresh).mean(),
                 s["mean"], s["mde_multiple"]))
        exp2_out[key] = row_out
    best_m_avg = min(d["m_ladder"], key=lambda m: np.mean([r["exp2"]["m%d" % m]["avg_m"] for r in rows]))
    best_m_clust = min(d["m_ladder"], key=lambda m: np.mean([r["exp2"]["m%d" % m]["clust_m"] for r in rows]))
    print("\n  optimal m (plain avg_m): %d     optimal m (clust_m, k=3 combined): %d" % (best_m_avg, best_m_clust))
    out["exp2"] = exp2_out
    out["exp2_optimal_m"] = {"avg_m": best_m_avg, "clust_m": best_m_clust}

    # ============================================================= EXPERIMENT 3
    print("\n--- EXPERIMENT 3: WEIGHTED AVERAGING, NESTED CV ---\n")
    exp3_out = {}
    families = {
        "score": [("score_c%s" % c, c) for c in d["score_c_grid"]],
        "rank": [("rank_p%s" % p, p) for p in d["rank_p_grid"]],
        "dist": [("dist_c%s" % c, c) for c in d["dist_c_grid"]],
    }
    F = sorted(set(folds.tolist()))
    for fam, grid in families.items():
        keys = [k for k, _ in grid]
        M = np.array([[r["exp3"][k] for k in keys] for r in rows], float)  # (n, ngrid)
        # in-sample oracle over the grid (leakage ceiling)
        insample_idx = int(np.argmin(M.mean(0)))
        insample_val = float(M[:, insample_idx].mean())
        insample_param = grid[insample_idx][1]
        # nested CV: leave-one-fold-out
        achieved = np.zeros(len(rows))
        chosen_params = {}
        for f in F:
            train = folds != f; test = folds == f
            tr_idx = int(np.argmin(M[train].mean(0)))
            achieved[test] = M[test, tr_idx]
            chosen_params[f] = grid[tr_idx][1]
        s = _stat(achieved - inc, folds, rng)
        s_insample = _stat(M[:, insample_idx] - inc, folds, rng)
        print("  %s-weighted:" % fam)
        print("    in-sample-oracle param=%s  mean=%.4f  vs avg_75: %s   <- LEAKAGE CEILING, not achievable"
              % (insample_param, insample_val, _fmt(s_insample)))
        print("    nested-CV achieved     chosen params per fold=%s  mean=%.4f  vs avg_75: %s"
              % (chosen_params, achieved.mean(), _fmt(s)))
        flag = " <<< BEATS INCUMBENT (held-out)" if s["ci_iid"][1] < 0 and s["mde_multiple"] > 1.0 else ""
        if flag:
            print("   " + flag)
        exp3_out[fam] = dict(insample_param=insample_param, insample_mean=insample_val,
                              insample_vs_incumbent=s_insample, chosen_params_per_fold=chosen_params,
                              achieved_mean=float(achieved.mean()), achieved_vs_incumbent=s)
        print()
    out["exp3"] = exp3_out

    path = os.path.join(RES, "agentA_report.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)
    print("\nwrote", path)
    return out


if __name__ == "__main__":
    report()
