"""s22/results/d_rg_full126.py -- extend L27/L28's compactness-disagreement channel (rg_z,
rg_gap) from the n<=13 subset (75 targets, latent-restricted) to the FULL n=126 pool, and
attempt the Angstrom conversion BRIEF s22 names as the campaign's central open question:
does rg_z/rg_gap let a NATIVE-FREE, HELD-OUT-FOLD router beat the fixed incumbent (3.048)?
"""
import json, os, sys
import numpy as np

sys.path.insert(0, ".")
from s12 import instrument as I

CA_BOND = 3.8046

def rg_of(P):
    P = np.asarray(P, float)
    c = P - P.mean(0, keepdims=True)
    return float(np.sqrt((c*c).sum() / len(P)))

def run():
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        dhat = np.asarray(dg["expected"], float)
        d2 = float((dhat**2).sum() + (n-1)*CA_BOND**2)
        rg_disto = float(np.sqrt(d2 / (n*n)))
        rgs = np.array([rg_of(w) for w in W], float)
        rows.append({
            "pdb": pdb, "n": n, "fold": fold,
            "rg_disto": rg_disto,
            "rg_pool_mean": float(rgs.mean()),
            "rg_pool_sd": float(rgs.std(ddof=1)),
            "rg_gap": float(rg_disto - rgs.mean()),
            "rg_z": float((rg_disto - rgs.mean()) / max(rgs.std(ddof=1), 1e-9)),
        })
        if (c+1) % 20 == 0:
            print(f"  {c+1}/{len(tg)}", flush=True)
    json.dump({"rows": rows, "n_expected": len(tg), "complete": len(rows)==len(tg)},
              open("s22/results/d_rg_full126.json", "w"))
    return rows

if __name__ == "__main__":
    run()
