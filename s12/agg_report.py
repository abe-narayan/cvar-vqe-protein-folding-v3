"""Reporting helpers: ladder tables, paired stats, FAIL18 split, per-fold, concentration."""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I

TG = I.targets()
PDBS = [t["pdb"] for t in TG]
FOLDS = np.array([t["fold"] for t in TG])
ISF = np.array([p in set(I.FAIL18) for p in PDBS])


def vec(rows, arm):
    return np.array([rows[p][arm] for p in PDBS], float)


def line(name, x, base=None):
    x = np.asarray(x, float)
    o = {"arm": name, "mean": round(float(x.mean()), 4), "median": round(float(np.median(x)), 4),
         "FAIL18": round(float(x[ISF].mean()), 4), "other108": round(float(x[~ISF].mean()), 4),
         "frac<2": round(float((x < 2.0).mean()), 3)}
    if base is not None:
        st = I.paired(x, np.asarray(base, float), folds=FOLDS, names=PDBS)
        o["d"] = round(st["mean_diff"], 4)
        o["ci"] = [round(st["ci95"][0], 4), round(st["ci95"][1], 4)]
        o["W/L"] = f"{st['n_better']}/{st['n_worse']}"
        o["drop10"] = round(st["drop_top10_mean_diff"], 4)
        o["drop20"] = round(st["drop_top20_mean_diff"], 4)
        o["per_fold"] = {k: round(v, 3) for k, v in st["per_fold"].items()}
    return o


def table(rows, cols=None):
    hdr = ["arm", "mean", "median", "FAIL18", "other108", "frac<2", "d", "ci", "W/L", "drop10", "drop20"]
    out = []
    keys = cols or hdr
    out.append("| " + " | ".join(keys) + " |")
    out.append("|" + "---|" * len(keys))
    for r in rows:
        out.append("| " + " | ".join(str(r.get(k, "")) for k in keys) + " |")
    return "\n".join(out)
