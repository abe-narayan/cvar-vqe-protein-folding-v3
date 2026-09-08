"""Aggregation for obj_* result files: per-condition means, FAIL18 split, paired stats."""
from __future__ import annotations
import os, json, sys, collections
import numpy as np
from s12 import instrument as I
from s12 import obj_common as OC

RES = os.path.join(I.ROOT, "s12", "results")


def load(name):
    with open(os.path.join(RES, name if name.endswith(".json") else name + ".json")) as fh:
        return json.load(fh)


def frame(name, key="fit_rmsd"):
    """-> (order of pdbs, {tag: np.array aligned to order}, meta)."""
    obj = load(name)
    rows = obj["rows"]
    conds = obj["spec"]["conds"]
    tags = [c.get("tag", c["kind"]) for c in conds]
    pdbs = sorted({r["pdb"] for r in rows})
    idx = {p: k for k, p in enumerate(pdbs)}
    out = {}
    extra = collections.defaultdict(lambda: np.full(len(pdbs), np.nan))
    for r in rows:
        t = tags[r["cond"]]
        out.setdefault(t, np.full(len(pdbs), np.nan))[idx[r["pdb"]]] = r.get(key, np.nan)
        for q in ("mae", "r", "sub_best", "sub_mean", "argmin_rmsd", "avg_rmsd", "fit_rmsd", "lam_rmsd"):
            if q in r:
                extra[(t, q)][idx[r["pdb"]]] = r[q]
    return pdbs, out, dict(extra), conds


def table(name, key="fit_rmsd", ref=None):
    pdbs, out, extra, conds = frame(name, key)
    tgt = {t["pdb"]: t for t in OC.targets()}
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    folds = np.array([tgt[p]["fold"] for p in pdbs])
    lines = []
    for t, v in out.items():
        mae = np.nanmean(extra.get((t, "mae"), np.array([np.nan])))
        rr = np.nanmean(extra.get((t, "r"), np.array([np.nan])))
        row = dict(tag=t, mean=float(np.nanmean(v)), mae=float(mae), r=float(rr),
                   fail18=float(np.nanmean(v[f18])), other=float(np.nanmean(v[~f18])),
                   frac2=float(np.nanmean(v < 2.0)),
                   sub_best=float(np.nanmean(extra.get((t, "sub_best"), np.array([np.nan])))),
                   argmin=float(np.nanmean(extra.get((t, "argmin_rmsd"), np.array([np.nan])))))
        if ref is not None and ref in out:
            st = I.paired(v, out[ref], folds=folds, names=pdbs)
            row.update(d=st["mean_diff"], ci=st["ci95"], W=st["n_better"], L=st["n_worse"],
                       drop10=st["drop_top10_mean_diff"], per_fold=st["per_fold"])
        lines.append(row)
    return lines, pdbs, out, extra, folds, f18


def show(name, key="fit_rmsd", ref=None, sort=True):
    lines, *_ = table(name, key, ref)
    if sort:
        lines = sorted(lines, key=lambda r: r["mean"])
    hdr = f"{'tag':28s} {'MAE':>6s} {'r':>6s} {key:>8s} {'FAIL18':>7s} {'other':>7s} {'<2A':>5s} {'top75best':>9s} {'argmin':>7s}"
    print(hdr); print("-" * len(hdr))
    for r in lines:
        s = (f"{r['tag']:28s} {r['mae']:6.3f} {r['r']:6.3f} {r['mean']:8.3f} {r['fail18']:7.3f} "
             f"{r['other']:7.3f} {r['frac2']:5.2f} {r['sub_best']:9.3f} {r['argmin']:7.3f}")
        if "d" in r:
            s += f"   d {r['d']:+.3f} [{r['ci'][0]:+.3f},{r['ci'][1]:+.3f}] {r['W']}W/{r['L']}L drop10 {r['drop10']:+.3f}"
        print(s)
    return lines


if __name__ == "__main__":
    show(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "fit_rmsd",
         sys.argv[3] if len(sys.argv) > 3 else None)
