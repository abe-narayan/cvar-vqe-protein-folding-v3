"""LANE L -- the ladder verdict.  Consumes `L2_ladder_{long,short}.jsonl`, emits P1-P3.

PREREG `s32/PREREG_S32_L.md` @ 88f2da39, L-H2.  The three predictions were fixed before
any ladder number existed and are evaluated here mechanically, each to its own verdict:

    P1  pool headroom survives: (avg75 - pool_best) >= 0.75 A on `long40`
    P2  selection is still the largest of the deployable rungs
    P3  the projection cost still tracks how non-physical the projected object is:
        a REAL MEMBER projects for < 0.05 A, a dense 75-member average for more

Every arm is `stats_lib.compare` (LOWER IS BETTER) with the instrument's OWN folds.  The
MDE is printed beside every effect.  < 0.7x MDE is NOT A RESULT; 0.7-1.0x is NOT MEASURED.

    python -m s32.s32_L_analyse
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
RUNGS = ["pool_best", "sparse_s10", "top75_best", "avg75", "avg75_random"]


def _load(kind):
    p = os.path.join(RESULTS, f"L2_ladder_{kind}.jsonl")
    rows = []
    seen = set()
    for line in open(p):
        try:
            r = json.loads(line)
        except Exception:
            continue
        if r["pdb"] in seen:
            continue
        seen.add(r["pdb"]); rows.append(r)
    rows.sort(key=lambda r: r["pdb"])
    return rows


def _tab(rows, basis):
    return {k: np.array([r[basis][k] for r in rows], float) for k in RUNGS}


def main():
    from s24 import stats_lib as ST
    out = {"basis": "built chain (projector output), pre-AMBER", "arms": {}, "ladders": {}}
    for kind in ("short", "long"):
        p = os.path.join(RESULTS, f"L2_ladder_{kind}.jsonl")
        if not os.path.exists(p):
            print(f"!! {kind}: no rows yet"); continue
        rows = _load(kind)
        folds = np.array([r["fold"] for r in rows])
        names = [r["pdb"] for r in rows]
        ch = _tab(rows, "chain"); cl = _tab(rows, "cloud")
        lens = np.array([r["n"] for r in rows], float)
        # COMPLETENESS.  `core/project._outpath` records this failure twice: a run stopped
        # early wrote the canonical filename and a partial table over an easier subset
        # looked finished.  The expected count is stated here and carried in the artefact,
        # so a reader can tell a finished ladder from an interrupted one.
        want = (45 if kind == "long" else 126)
        lad = {"n": len(rows), "n_expected": want, "complete": len(rows) >= want,
               "mean_len": float(lens.mean()),
               "len_range": [int(lens.min()), int(lens.max())],
               "fold_sizes": {int(f): int((folds == f).sum()) for f in sorted(set(folds))}}
        for k in RUNGS:
            lad[k] = {"chain_mean": float(ch[k].mean()),
                      "chain_median": float(np.median(ch[k])),
                      "chain_se": float(ch[k].std(ddof=1) / len(rows) ** 0.5),
                      "cloud_mean": float(cl[k].mean()),
                      "projection_cost": float((ch[k] - cl[k]).mean())}
        # ---- the rungs, as differences, each with its own MDE
        cmp = {}
        for lab, a, b in [
            ("pool_headroom  avg75 - pool_best", "avg75", "pool_best"),
            ("sparse_gain    sparse_s10 - pool_best", "sparse_s10", "pool_best"),
            ("retrieval      top75_best - pool_best", "top75_best", "pool_best"),
            ("selection      avg75 - top75_best", "avg75", "top75_best"),
            ("filter_skill   avg75 - avg75_random", "avg75", "avg75_random"),
        ]:
            c = ST.compare(ch[a], ch[b], folds=folds, names=names, label=f"{kind}:{lab}")
            cmp[lab] = {kk: c[kk] for kk in
                        ("n", "mean_a", "mean_b", "effect", "median_effect", "se", "mde",
                         "effect_over_mde", "ci95_iid", "n_better", "n_worse")}
            if "per_fold" in c:
                cmp[lab]["per_fold"] = c["per_fold"]
            for kk in c:
                if kk.startswith("fold") and kk not in cmp[lab]:
                    cmp[lab][kk] = c[kk]
        lad["comparisons"] = cmp
        # ---- draw control discipline (contract rule 10)
        lad["draw_sd_chain_mean"] = float(np.mean([r["draw_sd"]["chain"] for r in rows]))
        out["ladders"][kind] = lad

    # ---- the registered predictions
    v = {}
    if "long" in out["ladders"]:
        L = out["ladders"]["long"]
        head = L["avg75"]["chain_mean"] - L["pool_best"]["chain_mean"]
        v["P1_pool_headroom"] = {
            "value": head, "threshold": 0.75,
            "verdict": "HOLDS" if head >= 0.75 else "FAILS",
            "meaning": ("generation is NOT the bottleneck at 40-60 residues"
                        if head >= 0.75 else
                        "generation IS the bottleneck at 40-60 residues")}
        rungs = {"retrieval (pool->top75)":
                 L["top75_best"]["chain_mean"] - L["pool_best"]["chain_mean"],
                 "selection (top75_best->avg75)":
                 L["avg75"]["chain_mean"] - L["top75_best"]["chain_mean"]}
        big = max(rungs, key=rungs.get)
        v["P2_selection_largest"] = {
            "rungs": rungs, "largest": big,
            "verdict": "HOLDS" if big.startswith("selection") else "FAILS"}
        pm = L["pool_best"]["projection_cost"]
        av = L["avg75"]["projection_cost"]
        v["P3_projection_cost"] = {
            "real_member_cost": pm, "dense_average_cost": av,
            "verdict": "HOLDS" if (abs(pm) < 0.05 and av > abs(pm)) else "FAILS",
            "note": "contract rule 16: the cost is a property of the object projected"}
    out["predictions"] = v
    out["complete"] = all(l["complete"] for l in out["ladders"].values())         and set(out["ladders"]) == {"short", "long"}
    if not out["complete"]:
        for d in v.values():
            d["verdict"] = "PARTIAL -- " + d["verdict"]
    path = os.path.join(RESULTS, "L2_ladder_verdict.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)

    # ---- print
    for kind in ("short", "long"):
        if kind not in out["ladders"]:
            continue
        L = out["ladders"][kind]
        tag = "" if L["complete"] else f"  ** PARTIAL {L['n']}/{L['n_expected']} **"
        print(f"\n===== {kind.upper()}  n={L['n']}  len {L['len_range'][0]}-"
              f"{L['len_range'][1]} (mean {L['mean_len']:.1f})   "
              f"BASIS: built chain, pre-AMBER ====={tag}")
        print(f"{'rung':<16}{'chain':>9}{'median':>9}{'SE':>8}{'cloud':>9}{'proj cost':>11}")
        for k in RUNGS:
            d = L[k]
            print(f"{k:<16}{d['chain_mean']:9.4f}{d['chain_median']:9.4f}"
                  f"{d['chain_se']:8.4f}{d['cloud_mean']:9.4f}"
                  f"{d['projection_cost']:+11.4f}")
        print("\n  differences (LOWER IS BETTER; effect = a - b):")
        for lab, c in L["comparisons"].items():
            flag = ("RESULT" if abs(c["effect_over_mde"]) >= 1.0 else
                    "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "NOT A RESULT")
            print(f"    {lab:<40} {c['effect']:+8.4f}  MDE {c['mde']:.4f}  "
                  f"{abs(c['effect_over_mde']):5.2f}x  W/L {c['n_better']}/{c['n_worse']}  {flag}")
    if v:
        print("\n===== REGISTERED PREDICTIONS (prereg 88f2da39) =====")
        for k, d in v.items():
            print(f"  {k}: {d['verdict']}")
            for kk, vv in d.items():
                if kk != "verdict":
                    print(f"      {kk}: {vv}")
    print("\nwrote", path)


if __name__ == "__main__":
    main()
