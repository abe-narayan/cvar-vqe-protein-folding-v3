"""FAIL18 forensics, step 3: matched controls, class enrichment, failure classification.

Builds a length+fold matched control group of 18 well-performing targets, recomputes
every dossier contrast against it (so length is not the explanation), runs Fisher exact
tests on the provenance flags, and classifies each of the 18 into the failure taxonomy
(a) query / (b) filter / (c) objective / (d) label / (e) other.
"""
from __future__ import annotations
import os, sys, json, math, itertools
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I

RES = os.path.join(ROOT, "s12", "results")


def load():
    doss = json.load(open(os.path.join(RES, "fail_dossier.json")))
    hdr = json.load(open(os.path.join(RES, "fail_headers.json")))["per_target"]
    for d in doss:
        d["hdr"] = hdr[d["pdb"]]
    return doss


def fisher(a, b, c, d):
    """2x2 exact test, two-sided by the sum-of-small-probabilities rule."""
    from math import comb
    n = a + b + c + d
    r1, r2, c1 = a + b, c + d, a + c
    def p(k):
        return comb(r1, k) * comb(r2, c1 - k) / comb(n, c1)
    p0 = p(a)
    lo = max(0, c1 - r2); hi = min(r1, c1)
    return float(sum(p(k) for k in range(lo, hi + 1) if p(k) <= p0 * (1 + 1e-9)))


def matched_controls(doss, k=18, seed=0):
    """For each FAIL18 target pick the best-performing unused other-108 target with the
    same length if possible (else nearest length), preferring the same fold."""
    fail = [d for d in doss if d["fail18"]]
    pool = [d for d in doss if not d["fail18"]]
    # 'works well': emitted rmsd in the good half
    pool = sorted(pool, key=lambda d: d["rmsd_emitted"])
    used, ctrl = set(), []
    for f in sorted(fail, key=lambda d: -d["n"]):
        cands = [c for c in pool if c["pdb"] not in used]
        # rank: exact-length & same fold, then exact length, then |dn|
        cands.sort(key=lambda c: (abs(c["n"] - f["n"]), c["fold"] != f["fold"],
                                  c["rmsd_emitted"]))
        best = cands[0]
        used.add(best["pdb"])
        ctrl.append(dict(fail=f["pdb"], ctrl=best["pdb"], n_f=f["n"], n_c=best["n"],
                         fold_f=f["fold"], fold_c=best["fold"],
                         rmsd_f=f["rmsd_emitted"], rmsd_c=best["rmsd_emitted"]))
    return ctrl


def stat(rows, path):
    vals = []
    for r in rows:
        v = r
        for p in path.split("."):
            v = v[p]
        vals.append(float(v))
    a = np.asarray(vals, float)
    return a


def boot_diff(a, b, n_boot=4000, seed=0):
    rng = np.random.default_rng(seed)
    d = float(np.nanmean(a) - np.nanmean(b))
    bs = np.array([np.nanmean(a[rng.integers(0, len(a), len(a))]) -
                   np.nanmean(b[rng.integers(0, len(b), len(b))]) for _ in range(n_boot)])
    return d, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def classify(d):
    """Failure taxonomy for one target (uses ORACLE labels: diagnostic only)."""
    tags = []
    # (a) query: the near-native band of the UNIVERSE is largely outside the K=500 pool
    if d["o_pool_best"] - d["o_uni_best"] > 0.5:
        tags.append("a_query")
    # (c) objective: the objective actively prefers the pool's bad members
    if d["o_native_score_pct"] > 0.5 or d["o_rho_score_rmsd"] < 0.25:
        tags.append("c_objective")
    # (b) filter (residual): band present in pool but excluded, without (c)
    if d["o_band_in_top75"] == 0 and "c_objective" not in tags:
        tags.append("b_filter")
    h = d["hdr"]
    if h["fibril"] or h["lasso"] or h["cyclic"] or h["bound"]:
        tags.append("d_label_context")
    return tags


def main():
    doss = load()
    by = {d["pdb"]: d for d in doss}
    fail = [d for d in doss if d["fail18"]]
    other = [d for d in doss if not d["fail18"]]

    ctrl_map = matched_controls(doss)
    ctrl = [by[c["ctrl"]] for c in ctrl_map]
    print("== matched controls (length/fold matched, best emitted RMSD) ==")
    for c in ctrl_map:
        print(f"  {c['fail']} n={c['n_f']:2d} f{c['fold_f']} rmsd={c['rmsd_f']:5.2f}   <->   "
              f"{c['ctrl']} n={c['n_c']:2d} f{c['fold_c']} rmsd={c['rmsd_c']:5.2f}")

    keys = ["n", "o_pool_best", "o_uni_best", "o_top75_best", "o_pool_mean",
            "o_band_size", "o_uband_in_pool", "o_band_score_pct",
            "o_native_score_pct", "o_rho_score_rmsd", "o_rho_score_rmsd_band",
            "o_argmin_rmsd", "rmsd_fit", "rmsd_emitted", "nat_rg", "nat_rg_over_ref",
            "ee_over_contour", "rg_gap_nat_pool",
            "distogram.mae", "distogram.bias", "distogram.r", "distogram.slope",
            "distogram.shells.s2_4.mae", "distogram.shells.s5_8.mae",
            "distogram.shells.s9p.mae",
            "distogram.shells.s2_4.bias", "distogram.shells.s5_8.bias",
            "distogram.shells.s9p.bias",
            "composition.frac_org_pool", "composition.frac_org_band",
            "composition.frac_org_top75", "composition.sim_pool_mean",
            "composition.n_distinct_seq_pool",
            "nat_ss_frac.H", "nat_ss_frac.E", "nat_ss_frac.C"]
    print(f"\n{'quantity':34s} {'FAIL18':>8s} {'MATCH18':>8s} {'other108':>8s} {'diff':>8s} {'ci95':>18s}")
    table = {}
    for k in keys:
        try:
            a = stat(fail, k); b = stat(ctrl, k); c = stat(other, k)
        except KeyError:
            continue
        dd, lo, hi = boot_diff(a, b)
        table[k] = dict(fail18=float(np.nanmean(a)), match18=float(np.nanmean(b)),
                        other108=float(np.nanmean(c)), diff=dd, ci=[lo, hi])
        print(f"{k:34s} {np.nanmean(a):8.3f} {np.nanmean(b):8.3f} {np.nanmean(c):8.3f} "
              f"{dd:8.3f} [{lo:7.3f},{hi:7.3f}]")

    # ---- provenance enrichment, Fisher exact -------------------------------------
    print("\n== provenance enrichment: FAIL18 vs other-108 (Fisher exact) ==")
    enr = {}
    flags = ["fibril", "lasso", "cyclic", "membrane", "bound", "is_xray", "is_fragment"]
    def get(d, f):
        return bool(d["hdr"][f])
    for f in flags + ["fibril_or_lasso", "constrained_topology"]:
        if f == "fibril_or_lasso":
            g = lambda d: d["hdr"]["fibril"] or d["hdr"]["lasso"]
        elif f == "constrained_topology":
            g = lambda d: (d["hdr"]["lasso"] or d["hdr"]["cyclic"]
                           or d["hdr"]["n_ssbond"] > 0 or d["hdr"]["n_link"] > 0)
        else:
            g = lambda d, f=f: get(d, f)
        a = sum(g(d) for d in fail); b = len(fail) - a
        c = sum(g(d) for d in other); e = len(other) - c
        pv = fisher(a, b, c, e)
        enr[f] = dict(fail=a, fail_n=len(fail), other=c, other_n=len(other), p=pv)
        print(f"  {f:22s} {a:2d}/18  {c:3d}/108   p={pv:.2e}")

    # who are the non-FAIL fibril/lasso targets, and how do they do?
    print("\n== fibril/lasso targets NOT in FAIL18 (the within-class control) ==")
    for d in other:
        if d["hdr"]["fibril"] or d["hdr"]["lasso"]:
            print(f"  {d['pdb']} n={d['n']:2d} emitted={d['rmsd_emitted']:5.2f} "
                  f"poolbest={d['o_pool_best']:5.2f} dgMAE={d['distogram']['mae']:5.2f} "
                  f"natpct={d['o_native_score_pct']:.2f} | {d['hdr']['title'][:50]}")

    # ---- subgroup: fibril/lasso vs the rest of FAIL18 ----------------------------
    fl = [d for d in fail if d["hdr"]["fibril"] or d["hdr"]["lasso"]]
    nfl = [d for d in fail if not (d["hdr"]["fibril"] or d["hdr"]["lasso"])]
    print(f"\n== inside FAIL18: fibril/lasso (n={len(fl)}) vs rest (n={len(nfl)}) ==")
    for k in ["rmsd_emitted", "o_pool_best", "distogram.mae", "distogram.r",
              "o_native_score_pct", "nat_ss_frac.H", "nat_ss_frac.E", "nat_rg_over_ref"]:
        print(f"  {k:26s} {np.nanmean(stat(fl,k)):7.3f} {np.nanmean(stat(nfl,k)):7.3f}")

    # ---- per-target classification -----------------------------------------------
    print("\n== failure classification (ORACLE diagnostic) ==")
    cls = {}
    for d in sorted(fail, key=lambda x: x["pdb"]):
        t = classify(d)
        cls[d["pdb"]] = t
        print(f"  {d['pdb']} n={d['n']:2d} emit={d['rmsd_emitted']:5.2f} "
              f"pool={d['o_pool_best']:5.2f} uni={d['o_uni_best']:5.2f} "
              f"dgMAE={d['distogram']['mae']:5.2f} rho={d['o_rho_score_rmsd']:+.2f} "
              f"natpct={d['o_native_score_pct']:.2f} band={d['o_band_size']:3d} :: {','.join(t)}")
    cnt = {}
    for t in cls.values():
        for x in t:
            cnt[x] = cnt.get(x, 0) + 1
    print("  class counts:", cnt)
    # same classifier on the matched controls (the null: does it fire on good targets?)
    cls_c = {d["pdb"]: classify(d) for d in ctrl}
    cntc = {}
    for t in cls_c.values():
        for x in t:
            cntc[x] = cntc.get(x, 0) + 1
    print("  class counts on MATCH18 control:", cntc)

    I.write("fail_contrast", dict(controls=ctrl_map, table=table, enrichment=enr,
                                  classification=cls, classification_control=cls_c,
                                  class_counts=cnt, class_counts_control=cntc))


if __name__ == "__main__":
    main()
