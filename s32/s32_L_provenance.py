"""LANE L -- the adversarial check on the long pool's headroom, run BEFORE quoting it.

The first `long40` ladder rows show ORACLE pool-best members at 0.89-1.21 A on 52-58
residue targets, against a 75-member average at 5.9-9.7 A.  A pool headroom that large is
exactly the shape a LEAKAGE artefact takes: if the best member is a window of a structural
homolog that the sequence filter failed to remove, the rung is measuring the PDB's
redundancy and not the retrieval's reach.

`long40`'s candidate filter drops any source chain sharing a **verbatim 9-mer** with the
target.  That is the statistic with a measured null (`s32_L_leaknull`: real 0.330, shuffled
0.000), but it is a SEQUENCE filter, and two chains can be structurally near-identical with
no 9-mer in common.  So this module opens the winners and asks what they actually are:

    for the ORACLE best member of each target's K=500 pool, report
      - the source chain and window start
      - full-length sequence identity, source chain vs target        (longer-normalised)
      - window-vs-target identity at matched length                  (the comparable unit)
      - the longest verbatim common substring
      - the member's ORACLE CA-RMSD

A headroom built on windows whose identity to the target is at background level is real
retrieval reach.  One built on 60-90% identical windows is redundancy, and the rung must be
withdrawn.  **The threshold is stated before the numbers are read: if the MEDIAN window
identity of the ORACLE best members exceeds 0.40, lane L withdraws the long pool-headroom
claim.**

    python -m s32.s32_L_provenance run
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

WITHDRAW_AT = 0.40


def _lcs_substr(a, b):
    """Longest common CONTIGUOUS substring length (the verbatim measure)."""
    best = 0
    for L in range(min(len(a), len(b)), best, -1):
        subs = {a[i:i + L] for i in range(len(a) - L + 1)}
        if any(b[i:i + L] in subs for i in range(len(b) - L + 1)):
            return L
    return 0


def run(verbose=True):
    from core import data as cdata
    from core import geometry as geo
    from s12 import instrument as I
    from s32 import s32_L_ladder as LD
    from s32 import s32_L_corpus as cp

    corpus = cp.load()
    by_pdb = {p: a for a, p in enumerate(corpus["pdb"])}
    tg = json.load(open(os.path.join(RESULTS, "long40_manifest.json")))["targets"]
    rows = []
    for k, t in enumerate(tg):
        W, sim, src = LD._bank_long(t)
        p = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
        if not os.path.exists(p):
            p = glob.glob(os.path.join(BASE, "prots", t["pdb"].lower() + ".pdb"))[0]
        _, coords, _, _ = geo.native_coords_from_pdb(p)
        nat = np.asarray(coords["CA"], float)
        rr = I.kabsch_rmsd_batch(W, nat)
        a = int(np.argmin(rr))
        chain, start = src[a].rsplit("_", 1)
        cseq = corpus["seq"][by_pdb[chain]]
        wseq = cseq[int(start):int(start) + t["n"]]
        rows.append(dict(
            pdb=t["pdb"], n=int(t["n"]), rr=float(rr[a]), source=src[a],
            chain_identity=float(cdata.identity_many(t["seq"], [cseq])[0]),
            window_identity=float(cdata.identity_many(t["seq"], [wseq])[0]),
            longest_verbatim=int(_lcs_substr(t["seq"], wseq)),
            blosum_rank=int(a)))
        if verbose and (k + 1) % 10 == 0:
            print(f"  {k+1}/{len(tg)}", flush=True)

    wi = np.array([r["window_identity"] for r in rows])
    ci = np.array([r["chain_identity"] for r in rows])
    lv = np.array([r["longest_verbatim"] for r in rows])
    rr = np.array([r["rr"] for r in rows])
    out = dict(n=len(rows), withdraw_threshold=WITHDRAW_AT,
               window_identity=dict(mean=float(wi.mean()), median=float(np.median(wi)),
                                    p90=float(np.percentile(wi, 90)), max=float(wi.max())),
               chain_identity=dict(mean=float(ci.mean()), median=float(np.median(ci)),
                                   max=float(ci.max())),
               longest_verbatim=dict(mean=float(lv.mean()), median=float(np.median(lv)),
                                     max=int(lv.max())),
               oracle_pool_best_cloud=dict(mean=float(rr.mean()),
                                           median=float(np.median(rr))),
               verdict=("WITHDRAW: the long pool headroom is redundancy"
                        if float(np.median(wi)) > WITHDRAW_AT else
                        "HOLDS: the ORACLE best members are not sequence homologs"),
               rows=rows)
    path = os.path.join(RESULTS, "L2_provenance.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
        print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    run()
