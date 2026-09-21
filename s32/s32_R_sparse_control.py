#!/usr/bin/env python
"""S32 LANE R -- the control that separates "SPARSE is cheap to project" from
"ALIGNED is cheap to project".

EXPLORATORY.  Registered in s32/MULTIPLICITY.md as R-12 when it was written, after the
ladder null (R1) returned cos_align +0.40 on the ORACLE sparse rungs and -0.05 on
production.  It is not in PREREG_S32_R.md; it is the control that result demands.

THE QUESTION.  S32 contract rule 16 reads the O-ladder as "the projection cost is a
property of the OBJECT", with sparse convex combinations paying ~0 and dense averages
paying ~0.16, and the natural reading is that SPARSITY is what makes projection cheap.
But the ladder's sparse rungs are ORACLE: their weights were fitted to the native.  My R1
measurement says the two regimes sit at the SAME manifold distance (0.705 vs 0.844) and
differ only in the ALIGNMENT of the projection displacement with the native error.

So the two hypotheses make opposite predictions for a sparse combination whose weights the
native never touched:

    "SPARSITY is cheap"   -> a random s=10 average projects for ~0 too, cos ~ +0.4
    "ALIGNMENT is cheap"  -> a random s=10 average pays the FULL orthogonal tax,
                             sqrt(e^2+d^2)-e, with cos ~ 0, despite being just as sparse

CONTRACT RULE 6: the control is matched to the operator's own space -- same sparsity s=10,
same construction (superpose on the subset medoid, uniform mean), same projection call,
same job.  The ONLY thing that differs from the ladder's `sparse_pool_s10` is whether the
native chose the members.

ARMS, per target:
    RANDSPARSE_d{0,1,2}   10 members drawn uniformly from the K=500 pool, 3 pinned draws
                          (contract rule 10: a random control needs its own distribution)
    SCORESPARSE           the 10 highest-scoring members of the filtered pool, uniform mean.
                          DEPLOYABLE.  Reported as a MECHANISM DIAGNOSTIC of the projection
                          regime, NOT as an endpoint proposal: the `m` axis is closed by
                          S31 (leave-fold-out m +0.0075, 63W/63L) and this lane does not
                          reopen it.
    PROD                  the 75-member production average, same job, as the paired base.
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
os.makedirs(RESULTS, exist_ok=True)

S = 10
DRAWS = 3
SEED0 = 32_120_007            # pinned here, in the committed source


def avg_of(W):
    """The production averaging operator, restricted to a subset: superpose the subset on
    its own medoid and take the uniform coordinate mean.  Identical construction to
    `core.pipeline.average`, which is what makes the contrast a contrast."""
    P = I.pairwise_rmsd(W)
    b = int(np.argmin(P.sum(1) / max(P.shape[0] - 1, 1)))
    return I.superpose_batch(W, W[b]).mean(0)


def scored(C, ca, nat):
    e = float(I.ca_rmsd(C, nat))
    d = float(I.ca_rmsd(ca, C))
    ch = float(I.ca_rmsd(ca, nat))
    return dict(cloud_rmsd=e, d_to_cloud=d, chain_rmsd=ch,
                price_observed=ch - e,
                price_isotropic_null=float(np.hypot(e, d)) - e,
                cos_align=float((e ** 2 + d ** 2 - ch ** 2) / max(2.0 * e * d, 1e-12)),
                cloud_vbond_mean=float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean()),
                chain_vbond_mean=float(np.linalg.norm(np.diff(ca, axis=0), axis=1).mean()),
                cloud_rg=float(np.sqrt(((C - C.mean(0)) ** 2).sum(-1).mean())))


def one(pdb, verbose=True):
    t0 = time.time()
    with np.load(os.path.join(S29_STRUCTS, f"{pdb}.npz"), allow_pickle=True) as z:
        Cprod = np.asarray(z["prod"], float)
        seq, fold, n = str(z["seq"]), int(z["fold"]), int(z["n"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)                     # ORACLE: scoring only
    pool = I.pool_idx(u)
    Wpool = np.asarray(u["W"], float)[pool].astype(np.float32).astype(float)
    sub = np.asarray(I.shipped_record(pdb)["sub"], int)

    out = {"pdb": pdb, "n": n, "fold": fold, "s": S, "arms": {}}
    pr = I.project(Cprod, seq, fold)
    out["arms"]["PROD"] = scored(Cprod, np.asarray(pr["ca"], float), nat)

    for dd in range(DRAWS):
        rng = np.random.default_rng(SEED0 + 7919 * dd + (abs(hash(pdb)) % 100_000))
        pick = rng.choice(len(Wpool), size=min(S, len(Wpool)), replace=False)
        C = avg_of(Wpool[pick])
        p = I.project(C, seq, fold)
        out["arms"][f"RANDSPARSE_d{dd}"] = scored(C, np.asarray(p["ca"], float), nat)

    #: the 10 best-scoring members of the production filtered prefix -- deployable, and a
    #: REGIME DIAGNOSTIC only (the m axis is closed by S31 and is not reopened here)
    C = avg_of(Wpool[sub[:S]])
    p = I.project(C, seq, fold)
    out["arms"]["SCORESPARSE"] = scored(C, np.asarray(p["ca"], float), nat)

    out["secs"] = round(time.time() - t0, 2)
    if verbose:
        print("%-6s " % pdb + "  ".join(
            "%s d=%.3f price%+.4f cos%+.3f" % (k[:11], v["d_to_cloud"],
                                               v["price_observed"], v["cos_align"])
            for k, v in out["arms"].items()) + "  %.0fs" % out["secs"], flush=True)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    a = ap.parse_args(argv)
    pdbs = [t["pdb"] for t in I.targets()]
    mine = [p for k, p in enumerate(pdbs) if k % a.nshard == a.shard]
    path = os.path.join(RESULTS, f"s32_R_sparsectl_shard{a.shard}.jsonl")
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
        r = one(p)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
    print("shard %d done (%d targets)" % (a.shard, len(mine)), flush=True)


if __name__ == "__main__":
    main()
