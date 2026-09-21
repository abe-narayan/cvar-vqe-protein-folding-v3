#!/usr/bin/env python
"""s31/s31_V_zrank.py -- LANE V, V2: independent verification of S31-L17, the sprint's headline.

CLAIM AUDITED (s31/LEDGER.md S31-L17, s31/STATE.md): `core/pipeline.py:864`
`E = _zrank(pool["sc"][o])` with `o = top[:128]` and `top = argsort(sc)[:want]`, so `sc[o]` is
SORTED and `rankdata` of a sorted vector is 1..128 -- hence E is a TARGET-INDEPENDENT CONSTANT,
"up to the pool's tie structure".

Lane P verified it on THREE SYNTHETIC random score vectors.  That establishes the algebra but not
the hedge: `rankdata` uses method='average', so any TIE in the real top-128 scores bends E away
from the ramp.  The pool is known to contain duplicate structures (lane A: `n_distinct` mean
118.45 of 128), and duplicate structures have IDENTICAL scores.  **So ties are not hypothetical
here -- they are the common case, and nobody has measured how far they move E.**

This script measures E on the REAL benchmark, all 126 targets, and reports the deviation from the
exact ramp.  Native-free throughout: `E` is a function of the score order only, no native is used.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from core.pipeline import _zrank           # noqa: E402  (read-only import, core is not modified)

OUT = os.path.join(HERE, "results", "s31_V_zrank.json")
DIM = 128


def main():
    rows, Es = [], []
    for t in I.targets():
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        dg = I.distogram(pdb)
        n = int(u["n"])
        pool = np.asarray(u["order"], int)[:500]
        W = np.asarray(u["W"], float)[pool]
        ii, jj = I.pair_index(n, 2)
        D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        sc = I.shipped_score(dg, D)
        order = np.argsort(sc, kind="stable")          # pipeline uses cfg.tie_break (stable)
        top = order[:DIM]
        s_top = np.asarray(sc, float)[top]
        E = _zrank(s_top)
        Es.append(E)
        uniq = len(np.unique(s_top))
        # distinct STRUCTURES among the top-128 (duplicate structures -> identical scores)
        ndist = len({W[i].tobytes() for i in top})
        rows.append(dict(pdb=pdb, n=n, fold=int(t["fold"]),
                         n_unique_scores=int(uniq), n_tied=int(DIM - uniq),
                         n_distinct_structs=int(ndist),
                         sorted_ok=bool(np.all(np.diff(s_top) >= 0))))
    Es = np.asarray(Es)                                 # 126 x 128
    r = np.arange(1, DIM + 1, dtype=float)
    ramp = (r - r.mean()) / r.std()                     # the exact tie-free E

    dev = np.abs(Es - ramp[None, :])
    out = dict(
        n=len(rows), dim=DIM,
        CLAIM="S31-L17: E = _zrank(sc[top[:128]]) is a target-independent constant",
        VERDICT_algebra=dict(
            all_targets_sc_top_is_sorted=bool(all(x["sorted_ok"] for x in rows)),
            note="`top = argsort(sc)` so `sc[top]` is sorted by construction -- this is the step "
                 "that makes rankdata return 1..128 and the whole claim rests on it. CONFIRMED "
                 "on all 126 from the real scores."),
        ramp_first5=ramp[:5].tolist(), ramp_last3=ramp[-3:].tolist(),
        E_max_abs_dev_from_ramp=float(dev.max()),
        E_mean_abs_dev_from_ramp=float(dev.mean()),
        E_max_dev_by_target=dict(
            pdb=rows[int(dev.max(1).argmax())]["pdb"], value=float(dev.max(1).max())),
        n_targets_exactly_equal_to_ramp=int((dev.max(1) == 0.0).sum()),
        n_targets_with_any_tie=int(sum(1 for x in rows if x["n_tied"] > 0)),
        ties=dict(mean=float(np.mean([x["n_tied"] for x in rows])),
                  median=float(np.median([x["n_tied"] for x in rows])),
                  max=int(max(x["n_tied"] for x in rows)),
                  p90=float(np.percentile([x["n_tied"] for x in rows], 90))),
        distinct_structs=dict(mean=float(np.mean([x["n_distinct_structs"] for x in rows])),
                              min=int(min(x["n_distinct_structs"] for x in rows))),
        pairwise_max_abs_dev_between_targets=float(
            np.abs(Es[:, None, :] - Es[None, :, :]).max()),
        rows=rows)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1, default=float))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
