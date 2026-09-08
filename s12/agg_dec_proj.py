"""Project the learned set decoder's weighted averages -> EMITTED CA-RMSD."""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_features as AF
from s12 import agg_decoder as D

OUTD = os.path.join(ROOT, "s12", "cache", "agg_decproj")
os.makedirs(OUTD, exist_ok=True)


def run(name, shard=0, nshard=1):
    z = np.load(os.path.join(D.OUT, f"{name}.npz"))
    tg = [t for k, t in enumerate(I.targets()) if k % nshard == shard]
    path = os.path.join(OUTD, f"{name}_{shard}_{nshard}.json")
    cur = json.load(open(path)) if os.path.exists(path) else {}
    for t in tg:
        if t["pdb"] in cur:
            continue
        o = AF.load(t["pdb"])
        w = np.asarray(z[t["pdb"]], float); w = w / w.sum()
        C = (o["Asup"].astype(float) * w[:, None, None]).sum(0)
        out = I.project(C, t["seq"], t["fold"])
        cur[t["pdb"]] = {"fit": float(I.ca_rmsd(out["fit_ca"], o["nat"].astype(float))),
                         "raw": float(I.ca_rmsd(C, o["nat"].astype(float)))}
        json.dump(cur, open(path, "w"))
        print(t["pdb"], flush=True)


def collect(name):
    out = {}
    for f in os.listdir(OUTD):
        if f.startswith(name + "_"):
            out.update(json.load(open(os.path.join(OUTD, f))))
    I.write(f"agg_decproj_{name}", out)
    print(name, "n", len(out), "fit", np.mean([v["fit"] for v in out.values()]))


if __name__ == "__main__":
    if sys.argv[1] == "collect":
        collect(sys.argv[2])
    else:
        run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
