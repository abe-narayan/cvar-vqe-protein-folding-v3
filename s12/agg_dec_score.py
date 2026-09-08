"""Double-scoring of every set-decoder arm through THREE terminals + the set-mean law.

Coordinator rule (adversarial agent's operator transfer law): a decoder that works by better
RANKING is nearly invisible through the m=75 coordinate average, so every arm must be scored
through argmin, the m=25 average and the m=75 average.  The retained-set MEAN is reported as
the primary cheap surrogate (their law: emitted_avg ~ 1.13 * set_mean - 0.96, R^2 0.89).

For a weighting decoder the three terminals are read off the SAME weight vector:
    argmin   = the argmax-weight candidate            (its own CA-RMSD)
    m=25     = uniform coordinate average of the 25 highest-weight candidates
    m=75     = the weighted average over all 75 (what the arm was trained to emit)
Baseline for all three: the shipped distogram score used the same way.
"""
from __future__ import annotations
import os, sys, json, glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_features as AF
from s12 import agg_decoder as D

PDBS = [t["pdb"] for t in I.targets()]
FOLDS = np.array([t["fold"] for t in I.targets()])
ISF = np.array([p in set(I.FAIL18) for p in PDBS])


def terminals(o, order, w=None):
    """order: candidate indices best-first.  Returns the three terminals + set stats."""
    A = o["Asup"].astype(float); nat = o["nat"].astype(float); rr = o["rr"].astype(float)
    out = {"argmin": float(rr[order[0]])}
    for m in (10, 25, 75):
        idx = order[:m]
        out[f"avg{m}"] = float(I.ca_rmsd(A[idx].mean(0), nat))
        out[f"setmean{m}"] = float(rr[idx].mean())
        out[f"setbest{m}"] = float(rr[idx].min())
    if w is not None:
        ww = np.asarray(w, float); ww = ww / ww.sum()
        out["wavg75"] = float(I.ca_rmsd((A * ww[:, None, None]).sum(0), nat))
        out["setmean_w"] = float((ww * rr).sum())
    return out


def run(out="agg_dec_terminals"):
    data = {p: AF.load(p) for p in PDBS}
    res = {}
    # ---- baseline: the shipped distogram score, same three terminals
    base = {}
    for p in PDBS:
        o = data[p]
        base[p] = terminals(o, np.argsort(o["sc"].astype(float), kind="stable"))
    res["SHIPPED_score"] = base
    # ---- uniform reference (no ranking at all)
    res["uniform"] = {p: terminals(data[p], np.arange(len(data[p]["rr"])),
                                   np.ones(len(data[p]["rr"]))) for p in PDBS}
    # ---- every saved decoder arm
    for f in sorted(glob.glob(os.path.join(D.OUT, "*.npz"))):
        nm = os.path.basename(f)[:-4]
        z = np.load(f)
        if not set(PDBS) <= set(z.files):
            continue
        arm = {}
        for p in PDBS:
            o = data[p]
            w = np.asarray(z[p], float)
            arm[p] = terminals(o, np.argsort(-w, kind="stable"), w)
        res[nm] = arm
    I.write(out, res)
    return res


def report(res=None):
    if res is None:
        res = json.load(open(os.path.join(I.RESULTS, "agg_dec_terminals.json")))
    keys = ["argmin", "avg10", "avg25", "avg75", "wavg75", "setmean25", "setmean75", "setbest75"]
    b = res["SHIPPED_score"]
    print("%-16s" % "arm", " ".join("%9s" % k for k in keys))
    for nm in ["SHIPPED_score", "uniform"] + sorted(k for k in res if k not in ("SHIPPED_score", "uniform")):
        a = res[nm]
        row = []
        for k in keys:
            v = [a[p][k] for p in PDBS if k in a[p]]
            row.append("%9.4f" % np.mean(v) if v else "%9s" % "-")
        print("%-16s" % nm, " ".join(row))
    print("\npaired vs SHIPPED_score, per terminal (d, ci, W/L):")
    for nm in sorted(k for k in res if k not in ("SHIPPED_score", "uniform")):
        a = res[nm]
        line = ["%-16s" % nm]
        for k in ("argmin", "avg25", "avg75", "wavg75"):
            if k not in a[PDBS[0]]:
                line.append("%-30s" % "-"); continue
            va = np.array([a[p][k] for p in PDBS])
            vb = np.array([b[p][k] for p in PDBS]) if k in b[PDBS[0]] else np.array([b[p]["avg75"] for p in PDBS])
            s = I.paired(va, vb, folds=FOLDS)
            line.append("%s %+0.3f[%+0.3f,%+0.3f]%d/%d" % (k, s["mean_diff"], s["ci95"][0], s["ci95"][1],
                                                           s["n_better"], s["n_worse"]))
        print("  ".join(line))


if __name__ == "__main__":
    r = run()
    report(json.load(open(os.path.join(I.RESULTS, "agg_dec_terminals.json"))))
