"""ARM 4 (coordinator request) + the alternative-objective screen, run on the CHEAP proxy
(top-75 coordinate average, r=0.99 with the projected structure) plus ranking diagnostics.

For every scorer: emitted avg-RMSD, top-75 best/mean, argmin RMSD, Spearman(score, true
RMSD) over the whole K=500 pool, Spearman inside the near-native band (pool_best+1.5), and
the percentile of the NATIVE structure under the same score (the score evaluated at the
native distance vector) -- all on tuning126.
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from s12 import instrument as I
from s12 import obj_common as OC
from s12 import obj_alt as A


def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    if ra.std() < 1e-9 or rb.std() < 1e-9:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def native_row(d, name, lfo):
    """The score a scorer assigns to the NATIVE structure (diagnostic)."""
    dd = dict(d)
    dd["D"] = np.vstack([d["dtrue"][None, :], d["D"]])
    s = A.SCORERS[name](dd, **lfo)
    return float(s[0]), s[1:]


def run(names, out="obj_arm4"):
    tg = OC.targets()
    z = np.load(os.path.join(OC.CACHE, "obj_lfo.npz"))
    rows = []
    for k, t in enumerate(tg):
        d = OC.load(t["pdb"])
        lfo = {"weights": z[f"{t['fold']}/wsh"], "recal": z[f"{t['fold']}/recal"],
               "Sig": z[f"{t['fold']}/Sig"]}
        rr = d["rr"]; band = rr <= rr.min() + I.BAND
        for nm in names:
            snat, sc = native_row(d, nm, lfo)
            e = OC.emit(d, sc, lam=None)
            pct = float((sc < snat).mean())
            rows.append(dict(pdb=t["pdb"], fold=t["fold"], tag=nm, avg_rmsd=e["avg_rmsd"],
                             argmin_rmsd=e["argmin_rmsd"], sub_best=e["sub_best"],
                             sub_mean=e["sub_mean"],
                             rho=spearman(sc, rr), rho_band=spearman(sc[band], rr[band]),
                             nat_pct=pct, n_band=int(band.sum()),
                             band_recall=float(np.isin(np.where(band)[0], np.argsort(sc)[:75]).mean()),
                             pool_best=float(rr.min())))
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tg)}", flush=True)
    I.write(out, {"rows": rows, "scorers": names})
    report(rows)
    return rows


def report(rows):
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        by[r["tag"]].append(r)
    base = {r["pdb"]: r["avg_rmsd"] for r in by["l1"]} if "l1" in by else None
    hdr = (f"{'scorer':22s} {'avg':>7s} {'FAIL18':>7s} {'oth108':>7s} {'top75b':>7s} {'argmin':>7s} "
           f"{'rho':>6s} {'rhoBnd':>7s} {'natPct':>7s} {'recall':>7s}")
    print(hdr); print("-" * len(hdr))
    for tag, v in sorted(by.items(), key=lambda kv: np.mean([x["avg_rmsd"] for x in kv[1]])):
        f = lambda q, sel=None: float(np.nanmean([x[q] for x in v if sel is None or sel(x)]))
        print(f"{tag:22s} {f('avg_rmsd'):7.3f} {f('avg_rmsd', lambda x: x['pdb'] in I.FAIL18):7.3f} "
              f"{f('avg_rmsd', lambda x: x['pdb'] not in I.FAIL18):7.3f} {f('sub_best'):7.3f} "
              f"{f('argmin_rmsd'):7.3f} {f('rho'):6.3f} {f('rho_band'):7.3f} {f('nat_pct'):7.3f} "
              f"{f('band_recall'):7.3f}")


if __name__ == "__main__":
    names = sys.argv[1].split(",") if len(sys.argv) > 1 else [
        "bayes", "l1", "sep_resid", "sep_only", "global_offset", "sep_mix", "shell_std",
        "scalefree", "sdw", "contact", "contact_l1", "gram3", "shellw", "recal", "maha",
        "oracle", "sep_resid_oracle"]
    run(names, sys.argv[2] if len(sys.argv) > 2 else "obj_arm4")
