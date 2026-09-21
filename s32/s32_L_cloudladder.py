"""LANE L -- the CLOUD ladder at full instrument size, both lengths, in ONE script.

WHY THIS EXISTS.  The built-chain ladder (`s32_L_ladder`) is the primary arm, but it is
projection-bound: 482 s per target at L ~ 55.  The CLOUD rungs need no projection, so the
full instrument can be measured now and the chain ladder confirms it later on a subset of
the same targets.  Both bases are reported; **the cloud is a diagnostic and a cloud gain is
not an endpoint gain** (contract rule 0).

CONTRACT RULE 7 -- a claim combining two lanes' numbers is audited by nobody.  Everything
below is re-derived HERE, in one pass, from coordinates: nothing is quoted from another
lane's artefact and nothing is quoted from another sprint.  In particular the canonical
`pool best 1.7108` recorded in `s12/instrument.py` is RECOMPUTED rather than cited, so the
two lengths are produced by identical code on the same day.

    python -m s32.s32_L_cloudladder run
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
K = 500
M = 75


def run(verbose=True):
    from core import geometry as geo
    from s12 import instrument as I
    from s24 import stats_lib as ST
    from s32 import s32_L_ladder as LD
    from s32 import s32_L_corpus as cp

    out = {"basis": "CA POINT CLOUD -- a diagnostic, NOT the built-chain endpoint",
           "K": K, "M": M}
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
                fold = int(t["fold"])
            else:
                W, _, _ = LD._bank_long(t)
                p = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
                if not os.path.exists(p):
                    p = glob.glob(os.path.join(BASE, "prots",
                                               t["pdb"].lower() + ".pdb"))[0]
                _, coords, _, _ = geo.native_coords_from_pdb(p)
                nat = np.asarray(coords["CA"], float)
                fold = int(t["fold"])
            rr = I.kabsch_rmsd_batch(W, nat)                            # ORACLE labels
            P = I.pairwise_rmsd(W[:M])
            C, _ = I.coordinate_average(W[:M], P)
            rows.append(dict(pdb=t["pdb"], n=int(len(nat)), fold=fold,
                             pool_best=float(rr.min()),
                             top75_best=float(rr[:M].min()),
                             pool_mean=float(rr.mean()),
                             top75_mean=float(rr[:M].mean()),
                             avg75=float(I.ca_rmsd(C, nat))))
            if verbose and (k + 1) % 25 == 0:
                print(f"  {kind} {k+1}/{len(tg)}", flush=True)
        folds = np.array([r["fold"] for r in rows])
        g = {q: np.array([r[q] for r in rows]) for q in
             ("pool_best", "top75_best", "pool_mean", "top75_mean", "avg75")}
        d = {q: dict(mean=float(v.mean()), median=float(np.median(v)),
                     se=float(v.std(ddof=1) / len(v) ** 0.5)) for q, v in g.items()}
        d["n"] = len(rows)
        d["mean_len"] = float(np.mean([r["n"] for r in rows]))
        cmp = {}
        for lab, a, b in [("headroom  avg75 - pool_best", "avg75", "pool_best"),
                          ("retrieval top75_best - pool_best", "top75_best", "pool_best"),
                          ("readout   avg75 - top75_best", "avg75", "top75_best")]:
            c = ST.compare(g[a], g[b], folds=folds, label=f"{kind} cloud: {lab}")
            cmp[lab] = {kk: c[kk] for kk in
                        ("n", "mean_a", "mean_b", "effect", "median_effect", "se", "mde",
                         "effect_over_mde", "n_better", "n_worse")}
            if "per_fold" in c:
                cmp[lab]["per_fold"] = c["per_fold"]
        d["comparisons"] = cmp
        d["rows"] = rows
        out[kind] = d
    path = os.path.join(RESULTS, "L2_cloud_ladder.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        for kind in ("short", "long"):
            d = out[kind]
            print(f"\n== {kind}  n={d['n']}  L~{d['mean_len']:.1f}   "
                  f"BASIS: CA CLOUD (diagnostic, not the endpoint)")
            for q in ("pool_best", "top75_best", "top75_mean", "pool_mean", "avg75"):
                print(f"   {q:<12}{d[q]['mean']:8.4f}  median {d[q]['median']:7.4f}"
                      f"  SE {d[q]['se']:.4f}")
            for lab, c in d["comparisons"].items():
                band = ("RESULT" if abs(c["effect_over_mde"]) >= 1 else
                        "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7
                        else "NOT A RESULT")
                print(f"   {lab:<34}{c['effect']:+8.4f}  MDE {c['mde']:.4f}  "
                      f"{abs(c['effect_over_mde']):5.2f}x  W/L {c['n_better']}/"
                      f"{c['n_worse']}  {band}")
        print("\nwrote", path, flush=True)
    return out


if __name__ == "__main__":
    run()
