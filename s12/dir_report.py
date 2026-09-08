"""Reporting for the magnitude ladder (dir_emitavg.json) -- the driver's own printer hit a
key collision after the data was already written, so the tables are rebuilt here."""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I

LADDER = (1.0, 2.0, 3.0, 4.0, 6.0)
ARMS = ("o_sign", "head", "head_soft", "rand_acc", "sep", "shellmaj", "const", "shuf",
        "o_tmaj", "head_tmaj")


def main():
    rows = json.load(open(os.path.join(ROOT, "s12", "results", "dir_emitavg.json")))
    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    out = {}
    for g, rs in grp.items():
        folds = np.array([r["fold"] for r in rs]); nm = [r["pdb"] for r in rs]
        pt = np.array([r["cells"]["pt"] for r in rs])
        by = np.array([r["cells"]["bayes"] for r in rs])
        tab = {"_pt": float(pt.mean()), "_bayes": float(by.mean()), "n": len(rs)}
        for D in LADDER:
            for a in ARMS:
                x = np.array([r["cells"][f"{a}@{D}"] for r in rs])
                s1 = I.paired(x, pt, folds=folds); s2 = I.paired(x, by, folds=folds, names=nm)
                tab[f"{a}@{D}"] = dict(
                    mean=float(x.mean()),
                    d_pt=s1["mean_diff"], ci_pt=s1["ci95"], wl_pt=[s1["n_better"], s1["n_worse"]],
                    d_bayes=s2["mean_diff"], ci_bayes=s2["ci95"],
                    wl_bayes=[s2["n_better"], s2["n_worse"]],
                    drop10=s2["drop_top10_mean_diff"], drop20=s2["drop_top20_mean_diff"],
                    per_fold=s2.get("per_fold"))
        out[g] = tab
    I.write("dir_emitavg_report", out)

    for g in grp:
        t = out[g]
        print(f"\n== {g} (n={t['n']}) == pt={t['_pt']:.3f} bayes={t['_bayes']:.3f}"
              f"   (delta emitted vs BAYES, coordinate average)")
        print(f"{'delta':>5s}" + "".join(f"{h:>11s}" for h in ARMS))
        for D in LADDER:
            print(f"{D:5.1f}" + "".join(f"{t[f'{a}@{D}']['d_bayes']:11.3f}" for a in ARMS))
        print(f"{'W/L':>5s}" + "".join(
            f"{t[f'{a}@4.0']['wl_bayes'][0]:5d}/{t[f'{a}@4.0']['wl_bayes'][1]:<5d}" for a in ARMS))
        print(f"{'d10':>5s}" + "".join(f"{t[f'{a}@4.0']['drop10']:11.3f}" for a in ARMS))
        print("  CI(head@4)  ", [round(v, 3) for v in t["head@4.0"]["ci_bayes"]],
              "  CI(rand_acc@4)", [round(v, 3) for v in t["rand_acc@4.0"]["ci_bayes"]],
              "  CI(o_sign@4)", [round(v, 3) for v in t["o_sign@4.0"]["ci_bayes"]])
        print("  per-fold head@4:", {k: round(v, 3) for k, v in t["head@4.0"]["per_fold"].items()})


if __name__ == "__main__":
    main()
