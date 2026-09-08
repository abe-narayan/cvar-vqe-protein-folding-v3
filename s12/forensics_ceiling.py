"""forensics: the per-group ceilings.  Where do the 18 land if each stage is made perfect?

ORACLE/DIAGNOSTIC throughout -- these are bounds, not methods.  Outcome `rmsd_avg` (see section 9
of the findings for the proxy fidelity).  Stages, each replacing exactly one component:
  base            production BLOSUM-500 -> shipped score top-75 -> coordinate average
  ofilt75         perfect FILTER: the 75 lowest-rr members of the shipped K=500 pool
  ofilt25         perfect FILTER at m=25
  okey_ofilt25    perfect KEY and FILTER: the 25 lowest-rr windows of the WHOLE universe
  obest           the single best window of the pool (the argmin ceiling already in the record)
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from s12 import instrument as I
from s12 import forensics_lib as L

if __name__ == "__main__":
    tg = I.targets(); rows = {}
    for k, t in enumerate(tg):
        pdb = t["pdb"]; u = I.load_univ(pdb); rr = u["rr"]; dg = I.distogram(pdb)
        pool = I.pool_idx(u); n = t["n"]; i, j = I.pair_index(n)
        D = I.pair_dists(u["W"][pool], i, j).astype(np.float32).astype(float)
        sc = I.shipped_score(dg, D)
        r = {}
        base = np.argsort(sc, kind="stable")[:I.M]
        r["base"] = I.ca_rmsd(I.coordinate_average(u["W"][pool[base]])[0], u["nat_ca"])
        for m, nm in ((75, "ofilt75"), (25, "ofilt25")):
            sel = pool[np.argsort(rr[pool], kind="stable")[:m]]
            r[nm] = I.ca_rmsd(I.coordinate_average(u["W"][sel])[0], u["nat_ca"])
        sel = np.argsort(rr, kind="stable")[:25]
        r["okey_ofilt25"] = I.ca_rmsd(I.coordinate_average(u["W"][sel])[0], u["nat_ca"])
        r["obest"] = float(rr[pool].min()); r["obest_univ"] = float(rr.min())
        rows[pdb] = r
        print(f"[{k+1:3d}/126] {pdb}", flush=True)
    names = [t["pdb"] for t in tg]; F = np.array([p in I.FAIL18 for p in names])
    agg = {}
    for m in ("base", "ofilt75", "ofilt25", "okey_ofilt25", "obest", "obest_univ"):
        v = np.array([rows[p][m] for p in names], float)
        agg[m] = dict(all=float(v.mean()), fail18=float(v[F].mean()), other=float(v[~F].mean()),
                      frac_under_2=float((v < 2.0).mean()))
    I.write("forensics_ceiling", dict(rows=rows, aggregate=agg))
    print(f"{'arm':14s} {'all':>8s} {'FAIL18':>8s} {'other108':>9s} {'<2A':>6s}")
    for m, v in agg.items():
        print(f"{m:14s} {v['all']:8.3f} {v['fail18']:8.3f} {v['other']:9.3f} {v['frac_under_2']:6.2f}")
