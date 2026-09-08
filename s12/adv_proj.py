"""ADVERSARIAL AUDIT 0c -- does instrument.project reproduce the production projection?"""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I

def main(k=8, seed=0):
    tg = I.targets()
    rng = np.random.default_rng(seed)
    pick = sorted(rng.choice(len(tg), k, replace=False).tolist())
    rows = []
    for ix in pick:
        t = tg[ix]; pdb = t["pdb"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        rec = I.shipped_record(pdb); sub = np.asarray(rec["sub"], int)
        C, b = I.coordinate_average(u["W"][p][sub])
        out = I.project(C, t["seq"], t["fold"])
        r = dict(pdb=pdb, n=t["n"], fold=t["fold"],
                 fit_maxabs=float(np.abs(out["fit_ca"] - np.asarray(rec["fit_ca"], float)).max()),
                 fit_rmsd=float(I.ca_rmsd(out["fit_ca"], np.asarray(rec["fit_ca"], float))),
                 arm_maxabs=float(np.abs(out["ca"] - np.asarray(rec["ca"], float)).max()),
                 arm_rmsd=float(I.ca_rmsd(out["ca"], np.asarray(rec["ca"], float))),
                 phi_maxabs=float(np.abs(out["phi"] - np.asarray(rec["phi"], float)).max()),
                 prod_rmsd_fit=float(rec["rmsd_fit"]),
                 mine_rmsd_fit=float(I.ca_rmsd(out["fit_ca"], u["nat_ca"])),
                 prod_rmsd_arm=float(rec["rmsd_arm"]),
                 mine_rmsd_arm=float(I.ca_rmsd(out["ca"], u["nat_ca"])))
        rows.append(r); print(json.dumps(r), flush=True)
    agg = dict(n=len(rows),
               max_fit_rmsd=float(max(x["fit_rmsd"] for x in rows)),
               max_arm_rmsd=float(max(x["arm_rmsd"] for x in rows)),
               max_dfit=float(max(abs(x["prod_rmsd_fit"]-x["mine_rmsd_fit"]) for x in rows)),
               max_darm=float(max(abs(x["prod_rmsd_arm"]-x["mine_rmsd_arm"]) for x in rows)))
    print(json.dumps(agg, indent=1)); I.write("adv_projection", dict(agg=agg, rows=rows))

if __name__ == "__main__":
    main()
