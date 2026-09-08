"""Report dir_emit.json (the PROJECTED confirmation) over however many targets are done."""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I


def main():
    rows = json.load(open(os.path.join(ROOT, "s12", "results", "dir_emit.json")))
    grp = {"done": rows, "fail18": [r for r in rows if r["fail18"]],
           "other": [r for r in rows if not r["fail18"]]}
    keys = [k for k in rows[0]["arms"] if not np.isnan(rows[0]["arms"][k])]
    out = {"n_done": len(rows), "delta": 4.0, "groups": {}}
    for g, rs in grp.items():
        if not rs:
            continue
        folds = np.array([r["fold"] for r in rs]); nm = [r["pdb"] for r in rs]
        tab = {}
        print(f"\n== {g} (n={len(rs)}) projected, delta=4.0 ==")
        print(f"{'arm':11s}{'avg':>8s}{'fit':>8s}{'lam':>8s}{'d_fit vs bayes':>16s}"
              f"{'CI':>20s}{'W/L':>10s}{'d10':>8s}{'MAE':>7s}{'r_sep':>7s}{'t75b':>7s}")
        for a in keys:
            fit = np.array([r["arms"][a] for r in rs], float)
            if np.isnan(fit).any():
                continue
            by = np.array([r["arms"]["bayes"] for r in rs], float)
            s = I.paired(fit, by, folds=folds, names=nm)
            row = dict(avg=float(np.mean([r["avg"][a] for r in rs])),
                       fit=float(fit.mean()),
                       lam=float(np.nanmean([r["lam"][a] for r in rs])),
                       d=s["mean_diff"], ci=s["ci95"], wl=[s["n_better"], s["n_worse"]],
                       drop10=s["drop_top10_mean_diff"] if s["drop_top10_mean_diff"] is not None else float("nan"), per_fold=s.get("per_fold"),
                       mae=float(np.mean([r["meta"][a]["mae"] for r in rs])),
                       r_sep=float(np.nanmean([r["meta"][a]["r_sep"] for r in rs])),
                       sub_best=float(np.mean([r["meta"][a]["sub_best"] for r in rs])))
            tab[a] = row
            print(f"{a:11s}{row['avg']:8.3f}{row['fit']:8.3f}{row['lam']:8.3f}{row['d']:16.3f}"
                  f"  [{row['ci'][0]:+.3f},{row['ci'][1]:+.3f}]"
                  f"{row['wl'][0]:5d}/{row['wl'][1]:<4d}{row['drop10']:8.3f}"
                  f"{row['mae']:7.3f}{row['r_sep']:7.3f}{row['sub_best']:7.3f}")
        out["groups"][g] = tab
    I.write("dir_emit_partial", out)


if __name__ == "__main__":
    main()
