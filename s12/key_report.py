"""Paired analysis of the retrieval-key arms produced by s12/key_emit.py."""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I

METRICS = ["fit", "proj", "argmin", "best", "mean", "band", "top_m_best", "top_m_mean",
           "band_recall_in_m"]


def load(names):
    rows = None
    for nm in names:
        p = os.path.join(I.RESULTS, nm + ".json")
        r = json.load(open(p))
        d = {x["pdb"]: x for x in r}
        if rows is None:
            rows = d
        else:
            for k, v in d.items():
                if k in rows:
                    rows[k]["arms"].update(v["arms"])
                else:
                    rows[k] = v
    return [rows[t["pdb"]] for t in I.targets() if t["pdb"] in rows]


def table(rows, arms, metric="fit"):
    out = {}
    for g, sel in [("all", lambda r: True), ("FAIL18", lambda r: r["f18"]),
                   ("other108", lambda r: not r["f18"])]:
        rs = [r for r in rows if sel(r) and all(a in r["arms"] for a in arms)]
        out[g] = {"n": len(rs)}
        for a in arms:
            out[g][a] = float(np.mean([r["arms"][a][metric] for r in rs]))
    return out


def compare(rows, a, b, metric="fit"):
    """paired a - b (a = candidate, b = incumbent). negative = candidate better."""
    res = {}
    for g, sel in [("all", lambda r: True), ("FAIL18", lambda r: r["f18"]),
                   ("other108", lambda r: not r["f18"])]:
        rs = [r for r in rows if sel(r) and a in r["arms"] and b in r["arms"]]
        if not rs:
            continue
        x = np.array([r["arms"][a][metric] for r in rs])
        y = np.array([r["arms"][b][metric] for r in rs])
        res[g] = I.paired(x, y, folds=np.array([r["fold"] for r in rs]),
                          names=[r["pdb"] for r in rs])
    return res


def fmt(cmp_, label):
    o = []
    for g in ("all", "FAIL18", "other108"):
        if g not in cmp_:
            continue
        c = cmp_[g]
        d10 = c["drop_top10_mean_diff"]; d20 = c["drop_top20_mean_diff"]
        o.append(f"{label:34s} {g:9s} n={c['n']:3d} {c['mean_a']:.3f} vs {c['mean_b']:.3f} "
                 f"d={c['mean_diff']:+.3f} CI[{c['ci95'][0]:+.3f},{c['ci95'][1]:+.3f}] "
                 f"W/L={c['n_better']}/{c['n_worse']} med={c['median_diff']:+.3f} "
                 f"d10={'  n/a ' if d10 is None else format(d10,'+.3f')} "
                 f"d20={'  n/a ' if d20 is None else format(d20,'+.3f')} "
                 f"folds={ {k: round(v,3) for k,v in c.get('per_fold',{}).items()} }")
    return "\n".join(o)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", required=True)
    ap.add_argument("--base", default="blosum")
    ap.add_argument("--metric", default="fit")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rows = load(a.files.split(","))
    arms = list(rows[0]["arms"])
    print(f"targets: {len(rows)}   arms: {arms}\n")
    for m in ("fit", "proj", "best", "mean", "band", "top_m_best", "argmin"):
        t = table(rows, arms, m)
        print(f"--- {m} ---")
        for g in ("all", "FAIL18", "other108"):
            print(f"  {g:9s} n={t[g]['n']:3d} " + "  ".join(f"{k}={t[g][k]:.3f}" for k in arms))
    print()
    allc = {}
    for arm in arms:
        if arm == a.base:
            continue
        for m in ("fit", "proj"):
            c = compare(rows, arm, a.base, m)
            allc[f"{arm}|{m}"] = c
            print(fmt(c, f"{arm} vs {a.base} [{m}]"))
        print()
    if a.out:
        I.write(a.out, {"table": {m: table(rows, arms, m) for m in
                                  ("fit", "proj", "best", "mean", "band", "top_m_best",
                                   "top_m_mean", "argmin", "band_recall_in_m")},
                        "compare": allc,
                        "per_target": {r["pdb"]: {k: v for k, v in r["arms"].items()} for r in rows}})
