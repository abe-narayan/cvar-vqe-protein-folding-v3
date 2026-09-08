"""Analysis of the ISO-MAE surface produced by `obj_run` on the corruption grid."""
from __future__ import annotations
import os, sys, json, collections
import numpy as np
from s12 import instrument as I
from s12 import obj_report as R


def main(name="obj_iso_cheap", key="avg_rmsd"):
    pdbs, out, extra, conds = R.frame(name, key)
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    # group by (family, mae), average over seeds
    grid = collections.defaultdict(list)
    ref = {}
    for tag, v in out.items():
        if "|" not in tag:
            ref[tag] = v
            continue
        fam, mae, seed = tag.split("|")
        grid[(fam, float(mae[3:]))].append(v)
    fams = sorted({k[0] for k in grid})
    maes = sorted({k[1] for k in grid})
    print(f"{'family':22s}" + "".join(f"{m:>9.3g}" for m in maes) + "   (emitted mean, avg-RMSD)")
    print("-" * (22 + 9 * len(maes)))
    rows = {}
    for fam in fams:
        vals = []
        for m in maes:
            g = grid.get((fam, m))
            vals.append(float(np.mean([x.mean() for x in g])) if g else np.nan)
        rows[fam] = vals
        print(f"{fam:22s}" + "".join(f"{v:9.3f}" for v in vals))
    for k, v in ref.items():
        print(f"{'REF ' + k:22s}" + f"{v.mean():9.3f}")
    # iso-MAE spread
    print("\nISO-MAE spread (max - min across families at fixed MAE):")
    for a, m in enumerate(maes):
        col = np.array([rows[f][a] for f in fams])
        good = ~np.isnan(col)
        best = fams[int(np.nanargmin(col))]; worst = fams[int(np.nanargmax(col))]
        print(f"  MAE {m:6.3g}:  {np.nanmin(col):.3f} ({best})  ..  {np.nanmax(col):.3f} ({worst})"
              f"   spread {np.nanmax(col) - np.nanmin(col):.3f}")
    # FAIL18 split at the real MAE
    print("\nAt the distogram's own MAE (2.339): FAIL18 / other-108")
    for fam in fams:
        g = grid.get((fam, 2.339))
        if not g:
            continue
        M = np.mean(g, 0)
        print(f"  {fam:22s} {M.mean():6.3f}   {M[f18].mean():6.3f} / {M[~f18].mean():6.3f}")
    for k, v in ref.items():
        print(f"  REF {k:18s} {v.mean():6.3f}   {v[f18].mean():6.3f} / {v[~f18].mean():6.3f}")
    I.write(name + "_summary", {"fams": fams, "maes": maes, "rows": rows,
                                "ref": {k: float(v.mean()) for k, v in ref.items()}})
    return rows, maes, fams, ref, out, extra, pdbs


if __name__ == "__main__":
    main(*(sys.argv[1:] or []))
