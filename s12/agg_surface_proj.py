"""Confirm headline cells of the operator x objective-quality surface through the REAL
L-BFGS projection (production stage 3b, lam=0 fit arm)."""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_common as A
from s12 import agg_surface as S

OUTD = os.path.join(ROOT, "s12", "cache", "agg_surfproj")
os.makedirs(OUTD, exist_ok=True)

CELLS = [("a0.0", 1), ("a0.0", 5), ("a0.0", 25), ("a0.0", 75),
         ("a0.5", 20), ("a0.5", 75), ("shipped", 5), ("shipped", 20), ("shipped", 75)]


def target_vec(d, o):
    pi, pj = d["pi"].astype(int), d["pj"].astype(int)
    dtrue = np.linalg.norm(d["nat"][pi] - d["nat"][pj], axis=-1)
    e = d["exp"].astype(float)
    a = float(o[1:])
    return np.maximum(dtrue + a * (e - dtrue), 1.5)


def run(shard, nshard):
    tg = [t for k, t in enumerate(I.targets()) if k % nshard == shard]
    path = os.path.join(OUTD, f"{shard}_{nshard}.json")
    cur = json.load(open(path)) if os.path.exists(path) else {}
    for t in tg:
        pdb = t["pdb"]
        rec = cur.get(pdb, {})
        d = A.load(pdb, mem=False)
        W = d["Wp"]; nat = d["nat"]
        pi, pj = d["pi"].astype(int), d["pj"].astype(int)
        D = I.pair_dists(W, pi, pj)
        P = S.pool_P(pdb, W)
        for o, m in CELLS:
            key = f"{o}|m{m}"
            if key in rec:
                continue
            sc = d["sc"].astype(float) if o == "shipped" else np.abs(D - target_vec(d, o)[None, :]).mean(1)
            idx = np.argsort(sc, kind="stable")[:m]
            C = S.subset_avg(W, P, idx)
            out = I.project(C, t["seq"], t["fold"])
            rec[key] = {"raw": float(I.ca_rmsd(C, nat)),
                        "fit": float(I.ca_rmsd(out["fit_ca"], nat)),
                        "lam": float(I.ca_rmsd(out["ca"], nat))}
        cur[pdb] = rec
        json.dump(cur, open(path, "w"))
        print(pdb, flush=True)


def collect():
    out = {}
    for f in os.listdir(OUTD):
        out.update(json.load(open(os.path.join(OUTD, f))))
    I.write("agg_surface_proj", out)
    print("n", len(out))


if __name__ == "__main__":
    if sys.argv[1] == "collect":
        collect()
    else:
        run(int(sys.argv[1]), int(sys.argv[2]))
