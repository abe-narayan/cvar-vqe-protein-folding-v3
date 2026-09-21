#!/usr/bin/env python
"""S32 LANE R -- does the isotropic-null FORMULA predict the S29 O-ladder's projection
prices from each rung's own manifold distance `d`?

PREREG: s32/PREREG_S32_R.md at 02754f5a (R1).  This is the registered R1 mechanism test
extended, at the coordinator's request, from one object to the whole ladder.

THE CLAIM UNDER TEST.  Kabsch RMSD is a metric on shape space, so for a cloud C with native
error e = RMSD(C, nat) and a projected chain X at manifold distance d = RMSD(X, C):

        if the displacement X - C is statistically ORTHOGONAL to the native error,
        RMSD(X, nat) = sqrt(e^2 + d^2),  i.e.  price = sqrt(e^2 + d^2) - e.

S32's contract rule 16 records five ladder rows whose prices range from -0.0052 to +0.1701
and calls the price "a property of the OBJECT".  If the formula holds, the object-dependence
is ENTIRELY through `d`, and one line predicts all five rows.  A control matched to the
operator's own displacement is contract rule 6; this is that control, applied per rung.

WHAT WOULD MAKE IT FAIL (contract rule 5): a rung whose observed price is materially BELOW
sqrt(e^2+d^2)-e means the projection moves TOWARD the native there (recoverable direction);
materially ABOVE means it moves away.  Both are visible per target and per rung.

ORACLE.  `best1_*`, `sparse_*` and `bestm*` are clouds S29 SELECTED USING THE NATIVE.  They
are used here only as GEOMETRIC OBJECTS to test a formula about projection; no endpoint claim
is made from them and every row carries the ORACLE flag.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
S29_STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")
S29_CHAIN = os.path.join(ROOT, "s29", "results", "s29_O_chain_rows.jsonl")
os.makedirs(RESULTS, exist_ok=True)

#: (npz key, the S29 ladder row it is, the price contract rule 16 records, is it ORACLE)
RUNGS = [
    ("prod",             "PRODUCTION 75-member uniform average",  +0.1622, False),
    ("bestm",            "dense prefix average, K=500 (m free)",  +0.1701, True),
    ("best1_pool",       "best single member, K=500",             -0.0030, True),
    ("sparse_pool_s10",  "sparse convex combination, K=500, s=10", +0.0002, True),
    ("sparse_pool_s20",  "sparse convex combination, K=500, s=20", -0.0052, True),
]


def s29_chain_rows():
    got = {}
    for f in [S29_CHAIN] + [os.path.join(os.path.dirname(S29_CHAIN), x)
                            for x in os.listdir(os.path.dirname(S29_CHAIN))
                            if x.startswith("s29_O_chain_rows") and x.endswith(".jsonl")]:
        if not os.path.exists(f):
            continue
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except Exception:                                        # noqa: BLE001
                    continue
                got[(r["pdb"], r["item"])] = r
    return got


def one(pdb, s29, verbose=True):
    t0 = time.time()
    with np.load(os.path.join(S29_STRUCTS, f"{pdb}.npz"), allow_pickle=True) as z:
        seq, fold, n = str(z["seq"]), int(z["fold"]), int(z["n"])
        clouds = {k: np.asarray(z[k], float) for k, _, _, _ in RUNGS if k in z.files}
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)                     # ORACLE: scoring only
    out = {"pdb": pdb, "n": n, "fold": fold, "rungs": {}}
    for key, name, _, oracle in RUNGS:
        C = clouds.get(key)
        if C is None or not np.isfinite(C).all():
            continue
        pr = I.project(C, seq, fold)
        ca = np.asarray(pr["ca"], float)
        e = float(I.ca_rmsd(C, nat))
        d = float(I.ca_rmsd(ca, C))
        obs = float(I.ca_rmsd(ca, nat))
        ref = s29.get((pdb, key), {}).get("rmsd_chain")
        out["rungs"][key] = dict(
            name=name, ORACLE_NOT_DEPLOYABLE=bool(oracle),
            cloud_rmsd=e, d_to_cloud=d, chain_rmsd=obs,
            price_observed=obs - e, price_isotropic_null=float(np.hypot(e, d)) - e,
            cloud_vbond_mean=float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean()),
            chain_vbond_mean=float(np.linalg.norm(np.diff(ca, axis=0), axis=1).mean()),
            s29_chain_rmsd=ref,
            s29_abs_diff=(abs(obs - ref) if ref is not None else None))
    out["secs"] = round(time.time() - t0, 2)
    if verbose:
        s = " ".join("%s d=%.3f obs%+.4f null%+.4f" %
                     (k[:9], v["d_to_cloud"], v["price_observed"], v["price_isotropic_null"])
                     for k, v in out["rungs"].items())
        print("%-6s %s  %.0fs" % (pdb, s, out["secs"]), flush=True)
    return out


def rows_path(shard):
    return os.path.join(RESULTS, f"s32_R_laddernull_shard{int(shard)}.jsonl")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    a = ap.parse_args(argv)
    s29 = s29_chain_rows()
    pdbs = [t["pdb"] for t in I.targets()]
    mine = [p for k, p in enumerate(pdbs) if k % a.nshard == a.shard]
    path = rows_path(a.shard)
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["pdb"])
                except Exception:                                        # noqa: BLE001
                    pass
    for p in mine:
        if p in done:
            continue
        r = one(p, s29)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
    print("shard %d done (%d targets)" % (a.shard, len(mine)), flush=True)


if __name__ == "__main__":
    main()
