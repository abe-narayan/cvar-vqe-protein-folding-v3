"""SPRINT 14 / VQE -- PART C analysis: the central crossing curve, with statistics.

Reads `s14/results/vqe_signal_<base>.json` and answers the sprint's question on both axes:

    at what objective quality does VQE/CVaR start to beat classical search,
    and does it ever beat it on STRUCTURE rather than on ENERGY?

Both crossings are reported separately because they need not coincide, and the difference
between them is the finding.

Every comparison goes through `s12.instrument.paired`: paired mean difference, bootstrap
95% CI, W/L, per-fold values, and drop-top-10 concentration. With nine enumerated targets a
drop-top-10 is vacuous, so drop-top-2 and drop-top-3 are reported instead and labelled as
such -- silently printing a None would be worse than saying the check cannot be run at this n.

    python -m s14.vqe_report [--base legacy|prior]
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from s12 import instrument as I
from s14 import vqe_lib as V

CLASSICAL = ("random", "greedy", "anneal", "ga")
QUANTUM = ("vqe_a1.0", "vqe_a0.25", "vqe_a0.1", "vqe_a0.05")
ALL_ARMS = CLASSICAL + QUANTUM


def load(base="legacy"):
    p = os.path.join(V.RESULTS, f"vqe_signal_{base}.json")
    d = json.load(open(p))
    rows = []
    for key, c in d["cells"].items():
        pdb, b, sig, bud = key.split("|")
        r = {"pdb": pdb, "signal": float(sig), "budget": int(bud),
             "rho": c["rho_obj_rmsd"], "space_best": c["space_best_rmsd"],
             "space_mean": c["space_mean_rmsd"],
             "exact_rmsd": c["exact"]["rmsd_returned"]}
        for a in ALL_ARMS:
            if a not in c:
                continue
            r[f"{a}|obj"] = c[a].get("objective_gap")
            r[f"{a}|rmsd"] = c[a].get("rmsd_returned")
            r[f"{a}|struct"] = c[a].get("structural_gap")
            r[f"{a}|best_seen"] = c[a].get("rmsd_best_seen")
            r[f"{a}|sel"] = (c[a].get("rmsd_returned", np.nan)
                             - c[a].get("rmsd_best_seen", np.nan))
            r[f"{a}|H"] = c[a].get("entropy_bits")
            r[f"{a}|div"] = c[a].get("diversity")
            r[f"{a}|p25"] = c[a].get("pmass_below_2.5")
            r[f"{a}|f25"] = c[a].get("frac_below_2.5")
        rows.append(r)
    return d, rows


def paired_small(a, b, names=None, folds=None):
    """`I.paired` plus drop-top-2/3, which is the honest concentration check at n=9."""
    o = I.paired(np.asarray(a), np.asarray(b), names=names, folds=folds)
    d = np.asarray(a, float) - np.asarray(b, float)
    order = np.argsort(d)
    o["drop_top2_mean_diff"] = float(d[order[2:]].mean()) if len(d) > 2 else None
    o["drop_top3_mean_diff"] = float(d[order[3:]].mean()) if len(d) > 3 else None
    o.pop("drop_top10_mean_diff", None)
    o.pop("drop_top20_mean_diff", None)
    return o


def table(rows, budget, field, title, fmt="{:8.3f}"):
    sigs = sorted({r["signal"] for r in rows})
    print(f"\n{title}  (budget = {budget} objective evaluations)")
    print(f"  {'signal':>7s} {'rho':>7s} {'exact':>7s} "
          + "".join(f"{a.replace('vqe_a','q'):>9s}" for a in ALL_ARMS))
    out = {}
    for s in sigs:
        rs = [r for r in rows if r["signal"] == s and r["budget"] == budget]
        if not rs:
            continue
        def col(a):
            v = [r.get(f"{a}|{field}") for r in rs]
            v = [np.nan if x is None else float(x) for x in v]
            return np.nan if all(np.isnan(v)) else float(np.nanmean(v))
        vals = [col(a) for a in ALL_ARMS]
        out[s] = dict(zip(ALL_ARMS, vals))
        print(f"  {s:7.2f} {np.mean([r['rho'] for r in rs]):+7.3f} "
              f"{np.mean([r['exact_rmsd'] for r in rs]):7.3f} "
              + "".join(fmt.format(v).rjust(9) for v in vals))
    return out


def crossing(rows, budget, base="legacy"):
    """The two crossings: objective axis and structural axis, with CIs."""
    sigs = sorted({r["signal"] for r in rows})
    print()
    print("=" * 104)
    print(f"THE CROSSING CURVE at budget {budget}: best VQE arm vs best classical arm")
    print("  paired over the nine enumerated targets; negative = VQE better")
    print("=" * 104)
    print(f"  {'signal':>7s} {'rho':>7s} | {'OBJECTIVE axis':>32s} | "
          f"{'STRUCTURE axis':>32s}")
    print(f"  {'':>7s} {'':>7s} | {'diff':>8s} {'CI95':>17s} {'W/L':>5s} | "
          f"{'diff':>8s} {'CI95':>17s} {'W/L':>5s}")
    out = {}
    for s in sigs:
        rs = sorted([r for r in rows if r["signal"] == s and r["budget"] == budget],
                    key=lambda r: r["pdb"])
        if len(rs) < 4:
            continue
        names = [r["pdb"] for r in rs]
        row = {}
        for axis, field in (("obj", "obj"), ("struct", "struct")):
            # best arm chosen by mean, per side -- predefined selection rule
            cq = min(QUANTUM, key=lambda a: np.nanmean([r[f"{a}|{field}"] for r in rs]))
            cc = min(CLASSICAL, key=lambda a: np.nanmean([r[f"{a}|{field}"] for r in rs]))
            p = paired_small([r[f"{cq}|{field}"] for r in rs],
                             [r[f"{cc}|{field}"] for r in rs], names=names)
            row[axis] = {"vqe_arm": cq, "classical_arm": cc, **p}
        out[s] = row
        o, t = row["obj"], row["struct"]
        print(f"  {s:7.2f} {np.mean([r['rho'] for r in rs]):+7.3f} | "
              f"{o['mean_diff']:+8.4f} [{o['ci95'][0]:+7.4f},{o['ci95'][1]:+7.4f}] "
              f"{o['n_better']:2d}/{o['n_worse']:<2d} | "
              f"{t['mean_diff']:+8.3f} [{t['ci95'][0]:+7.3f},{t['ci95'][1]:+7.3f}] "
              f"{t['n_better']:2d}/{t['n_worse']:<2d}")
    print("\n  A crossing is the first signal level at which the CI excludes zero on the")
    print("  favourable side.  The two axes are reported separately and never substituted.")
    return out


def budget_trap(rows, base="legacy"):
    """The known trap: is a VQE arm on the improving limb of a curve with a bad limit?"""
    print()
    print("=" * 104)
    print("THE BUDGET TRAP CHECK -- every arm against the CERTIFIED optimum")
    print("  'exact' is the certified global optimum of that exact objective over all")
    print("  262,144 configurations.  An arm that beats it has NOT solved the problem")
    print("  better; it has failed to optimise, and got lucky.")
    print("=" * 104)
    sigs = sorted({r["signal"] for r in rows})
    buds = sorted({r["budget"] for r in rows})
    print(f"  {'signal':>7s} {'exact':>7s} {'space best':>11s} {'random draw':>12s} "
          + "".join(f"{f'b={b}':>11s}" for b in buds) + "   (best arm's RMSD)")
    out = {}
    for s in sigs:
        line, ex, sb, sm = [], None, None, None
        for b in buds:
            rs = [r for r in rows if r["signal"] == s and r["budget"] == b]
            if not rs:
                line.append(np.nan); continue
            ex = np.mean([r["exact_rmsd"] for r in rs])
            sb = np.mean([r["space_best"] for r in rs])
            sm = np.mean([r["space_mean"] for r in rs])
            line.append(min(np.nanmean([r[f"{a}|rmsd"] for r in rs]) for a in ALL_ARMS))
        if ex is None:
            continue
        out[s] = {"exact": ex, "space_best": sb, "random": sm, "best_by_budget": line}
        print(f"  {s:7.2f} {ex:7.3f} {sb:11.3f} {sm:12.3f} "
              + "".join(f"{v:11.3f}" for v in line))
    print("\n  Where the budgeted columns are BETTER than 'exact', more search would make")
    print("  the answer worse -- the arm is on the improving limb of a curve whose limit")
    print("  is known to be bad.")
    return out


def main(argv):
    base = "legacy"
    if "--base" in argv:
        base = argv[argv.index("--base") + 1]
    d, rows = load(base)
    done = {(r["pdb"], r["signal"], r["budget"]) for r in rows}
    print(f"cells: {len(rows)}  targets: {len({r['pdb'] for r in rows})}"
          f"  base: {base}")
    out = {"config": d["config"], "n_cells": len(rows)}
    for bud in sorted({r["budget"] for r in rows}):
        print()
        print("=" * 104)
        print(f"BUDGET {bud}")
        print("=" * 104)
        out[f"obj_gap_{bud}"] = table(rows, bud, "obj",
                                      "OBJECTIVE GAP  E_found - E_exact (lower = better "
                                      "optimisation)", "{:8.4f}")
        out[f"struct_gap_{bud}"] = table(rows, bud, "struct",
                                         "STRUCTURAL GAP  RMSD_returned - RMSD_best_in_space")
        out[f"sel_gap_{bud}"] = table(rows, bud, "sel",
                                      "SELECTION GAP  RMSD_returned - RMSD_best_EVALUATED")
        out[f"pmass_{bud}"] = table(rows, bud, "p25",
                                    "PROBABILITY MASS below 2.5 A (VQE arms only "
                                    "have this)", "{:8.4f}")
        out[f"crossing_{bud}"] = crossing(rows, bud, base)
    out["budget_trap"] = budget_trap(rows, base)
    V.write(f"vqe_report_{base}", out)
    print(f"\nwritten -> s14/results/vqe_report_{base}.json")


if __name__ == "__main__":
    main(sys.argv[1:])
