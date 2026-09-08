"""s22/results/d_mreal_notrunc.py -- attack on mreal.py: does the transfer finding depend on the
DEGENERATE half-pool m=500(->250) rung?  Same seeding (SD.stable_rng(pdb,'s22mreal')) so the SAME
8 random half-splits are reproduced exactly, but the ladder used for SELECTION excludes m=500."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, ".")
from s12 import instrument as I
from s15 import seed as SD

MS_FULL = (500, 150, 75, 20, 5, 1)
MS_NOTRUNC = (150, 75, 20, 5, 1)   # drop the degenerate rung
REPS = 8
FIXED = 75

def run():
    tg = I.targets()
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        rng = SD.stable_rng(pdb, "s22mreal")   # IDENTICAL seed/stream to mreal.py
        rep = []
        for _r in range(REPS):
            perm = rng.permutation(len(W))     # consumes the RNG exactly as mreal.py does
            halves = (perm[:len(W)//2], perm[len(W)//2:])
            lad = []
            for h in halves:
                o = h[np.argsort(sc[h], kind="stable")]
                row = []
                for m in MS_NOTRUNC:
                    sel = o[:min(m, len(o))]
                    C, _b = I.coordinate_average(W[sel])
                    row.append(float(I.ca_rmsd(np.asarray(C, float), nat)))
                lad.append(row)
            A, B = np.array(lad[0]), np.array(lad[1])
            tiedA = np.flatnonzero(A <= A.min()+1e-12)
            tiedB = np.flatnonzero(B <= B.min()+1e-12)
            fixed_idx = MS_NOTRUNC.index(FIXED)
            rep.append({"selA_onB": float(B[tiedA].mean()), "selB_onA": float(A[tiedB].mean()),
                        "fixed_onB": float(B[fixed_idx]), "fixed_onA": float(A[fixed_idx])})
        agg = lambda k: float(np.mean([x[k] for x in rep]))
        rows.append({"pdb": pdb, "fold": int(t["fold"]),
                     "sel_heldout": 0.5*(agg("selA_onB")+agg("selB_onA")),
                     "fixed_heldout": 0.5*(agg("fixed_onB")+agg("fixed_onA"))})
        if (c+1) % 30 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"rows": rows, "ms": list(MS_NOTRUNC)}, open("s22/results/d_mreal_notrunc.json","w"))
    return rows

if __name__ == "__main__":
    rows = run()
    sel = np.array([r["sel_heldout"] for r in rows]); fix = np.array([r["fixed_heldout"] for r in rows])
    d = sel - fix
    se = d.std(ddof=1)/np.sqrt(len(d))
    print(f"\nWITHOUT the degenerate m=500(->250) rung, ladder = {MS_NOTRUNC}:")
    print(f"selected - fixed, held out: {d.mean():+.4f} SE {se:.4f} MDE {2.8016*se:.4f}  "
          f"W/L {int((d<0).sum())}/{int((d>0).sum())}")
