"""LANE L -- geometry secondaries for the one rung that CONSTRUCTS rather than selects.

Contract rule 15: an arm that constructs a structure rather than selecting one carries the
geometry checks from its first row.  Four of the five ladder rungs select a real deposited
window or a convex combination of posed windows; `avg75` CONSTRUCTS -- it is a uniform
coordinate average in the medoid frame, and a coordinate average is not a chain.

Project memory records that averaging CONTRACTS the backbone, and that the contraction is a
separation-dependent distortion rather than a uniform shrink
(`averaging-space-beats-the-objective`, where the originally quoted 25.8% was WITHDRAWN and
corrected to 3.5%).  That measurement was made at 9-16 residues.  Averaging 75 structures
that are mutually far apart must contract more, and at 40-60 residues they are much further
apart, so the cheap prediction is that the distortion grows with length.

Reported per target, on the PRE-PROJECTION cloud (the constructed object itself):

    rg_ratio        radius of gyration of the average / that of the native
    bond_mean       mean virtual CA-CA bond of the average (native backbone is ~3.80 A)
    bond_sd         its spread
    spread          mean pairwise CA-RMSD among the 75 members -- how far apart the
                    things being averaged actually are

    python -m s32.s32_L_geometry run
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
M = 75
K = 500


def _stats(C, nat, P):
    from core import geometry as geo
    C = np.asarray(C, float); nat = np.asarray(nat, float)
    b = np.linalg.norm(np.diff(C, axis=0), axis=1)
    bn = np.linalg.norm(np.diff(nat, axis=0), axis=1)
    iu = np.triu_indices(len(P), 1)
    return dict(rg=float(geo.radius_of_gyration(C)),
                rg_native=float(geo.radius_of_gyration(nat)),
                rg_ratio=float(geo.radius_of_gyration(C) /
                               geo.radius_of_gyration(nat)),
                bond_mean=float(b.mean()), bond_sd=float(b.std()),
                bond_native_mean=float(bn.mean()),
                bond_ratio=float(b.mean() / bn.mean()),
                spread=float(np.asarray(P)[iu].mean()))


def run(verbose=True):
    from core import geometry as geo
    from s12 import instrument as I
    from s32 import s32_L_ladder as LD
    from s32 import s32_L_corpus as cp

    out = {}
    for kind in ("short", "long"):
        rows = []
        if kind == "short":
            tg = I.targets()
        else:
            tg = json.load(open(os.path.join(RESULTS, "long40_manifest.json")))["targets"]
            cp.load()
        for k, t in enumerate(tg):
            if kind == "short":
                u = I.load_univ(t["pdb"])
                W = np.asarray(u["W"], float)[I.pool_idx(u, K)]
                nat = np.asarray(u["nat_ca"], float)
            else:
                W, _, _ = LD._bank_long(t)
                p = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
                if not os.path.exists(p):
                    p = glob.glob(os.path.join(BASE, "prots",
                                               t["pdb"].lower() + ".pdb"))[0]
                _, coords, _, _ = geo.native_coords_from_pdb(p)
                nat = np.asarray(coords["CA"], float)
            Ws = W[:M]
            P = I.pairwise_rmsd(Ws)
            C, _ = I.coordinate_average(Ws, P)
            r = _stats(C, nat, P)
            r.update(pdb=t["pdb"], n=int(len(nat)))
            rows.append(r)
            if verbose and (k + 1) % 25 == 0:
                print(f"  {kind} {k+1}/{len(tg)}", flush=True)
        agg = {}
        for key in ("rg_ratio", "bond_mean", "bond_sd", "bond_ratio", "spread",
                    "bond_native_mean"):
            v = np.array([r[key] for r in rows])
            agg[key] = dict(mean=float(v.mean()), median=float(np.median(v)),
                            se=float(v.std(ddof=1) / len(v) ** 0.5))
        agg["n"] = len(rows)
        agg["mean_len"] = float(np.mean([r["n"] for r in rows]))
        agg["rows"] = rows
        out[kind] = agg
    path = os.path.join(RESULTS, "L2_geometry_secondaries.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: {a: b for a, b in v.items() if a != "rows"}
                          for k, v in out.items()}, indent=1))
        print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    run()
