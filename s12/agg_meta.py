"""Adds pool-level metadata (org flag, BLOSUM sim, torsions) to the agg cache."""
from __future__ import annotations
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I

META = os.path.join(ROOT, "s12", "cache", "agg_meta")
os.makedirs(META, exist_ok=True)


def build():
    for t in I.targets():
        path = os.path.join(META, f"{t['pdb']}.npz")
        if os.path.exists(path):
            continue
        u = I.load_univ(t["pdb"])
        p = I.pool_idx(u)
        np.savez_compressed(path, org=np.asarray(u["org"])[p].astype(np.int8),
                            sim=np.asarray(u["sim"], np.float32)[p],
                            phi=np.asarray(u["PHI"], np.float32)[p],
                            psi=np.asarray(u["PSI"], np.float32)[p],
                            S=np.asarray(u["S"], np.int8)[p])
        del u
        print(t["pdb"], flush=True)


def load(pdb):
    z = np.load(os.path.join(META, f"{pdb}.npz"))
    return {k: z[k] for k in z.files}


if __name__ == "__main__":
    build()
    print("meta ok", len(os.listdir(META)))
