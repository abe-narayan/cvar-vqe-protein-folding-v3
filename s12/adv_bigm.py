"""Close the one gap in the m ladder: m = 300 / 500 through the REAL projection.

The coordinate-average arm says that on the FAIL18, discarding the distogram filter
entirely (averaging the whole K=500 BLOSUM pool) beats the shipped top-75 by -0.612 A.
That cell was never projected.  Run it on the FAIL18 and on a length/fold-matched control
of 18 non-FAIL18 targets so the contrast is not confounded by group composition.
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I

MS = [75, 150, 300, 500]

def main():
    tg = I.targets()
    fl = [t for t in tg if t["pdb"] in I.FAIL18]
    rest = [t for t in tg if t["pdb"] not in I.FAIL18]
    # length+fold matched control, greedy, one per FAIL18 target
    used, ctrl = set(), []
    for t in fl:
        cand = sorted([c for c in rest if c["pdb"] not in used],
                      key=lambda c: (abs(c["n"] - t["n"]), c["fold"] != t["fold"], c["pdb"]))
        ctrl.append(cand[0]); used.add(cand[0]["pdb"])
    rows = {}
    t0 = time.time()
    for grp, ts in (("FAIL18", fl), ("MATCH18", ctrl)):
        for t in ts:
            u = I.load_univ(t["pdb"]); p = I.pool_idx(u)
            W = np.asarray(u["W"][p], float); nat = u["nat_ca"]
            dg = I.distogram(t["pdb"], t["seq"], t["fold"])
            i, j = I.pair_index(t["n"])
            order = np.argsort(I.shipped_score(dg, I.pair_dists(W, i, j)), kind="stable")
            r = dict(grp=grp, pdb=t["pdb"], n=t["n"], fold=t["fold"])
            for m in MS:
                C, _ = I.coordinate_average(W[order[:m]])
                r[f"avg_m{m}"] = float(I.ca_rmsd(C, nat))
                r[f"fit_m{m}"] = float(I.ca_rmsd(I.project(C, t["seq"], t["fold"])["fit_ca"], nat))
            rows[t["pdb"]] = r
            print(grp, t["pdb"], " ".join(f"m{m}={r[f'fit_m{m}']:.2f}" for m in MS),
                  f"{time.time()-t0:.0f}s", flush=True)
    out = {"rows": rows}
    for grp in ("FAIL18", "MATCH18"):
        g = [r for r in rows.values() if r["grp"] == grp]
        b = np.array([r["fit_m75"] for r in g])
        out[grp] = {"n": len(g), "m75": float(b.mean())}
        for m in MS:
            if m == 75: continue
            a = np.array([r[f"fit_m{m}"] for r in g])
            s = I.paired(a, b, names=[r["pdb"] for r in g])
            out[grp][f"m{m}"] = float(a.mean())
            out[grp][f"paired_m{m}"] = s
            print(f"{grp} m{m}: {a.mean():.4f} vs m75 {b.mean():.4f}  d={s['mean_diff']:+.4f} "
                  f"CI[{s['ci95'][0]:+.4f},{s['ci95'][1]:+.4f}] W/L {s['n_better']}/{s['n_worse']}")
    I.write("adv_bigm", out)

if __name__ == "__main__":
    main()
