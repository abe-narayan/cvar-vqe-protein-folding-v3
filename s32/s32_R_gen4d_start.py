#!/usr/bin/env python
"""S32 LANE R -- EXPLORATORY R-16: the one branch signal that TRANSFERS, and why it buys
nothing.

Logged as R-16 in `s32/MULTIPLICITY.md`.  Not in `PREREG_S32_R.md`: it was written after the
best-of-K pricing showed that of six branch families only `GEN4D` has a split-half transfer
materially above zero.

WHY GEN4D IS DIFFERENT FROM THE OTHER FAMILIES.  `split_half_transfer` chooses a COLUMN on
half the targets and scores it on the other half, so a column has to mean the same thing on
every target.  For `GEN4D` a column is one NAMED generic start -- extended, alpha-helix,
beta-strand, polyproline II -- taken directly to lam = 0.3.  For `MEM75` a column is "the
r-th distogram-ranked member", and for `RAND0` the column index is arbitrary BY CONSTRUCTION,
which is exactly why `RAND0` is the zero-signal control for the transfer statistic.

THE DEPLOYABLE FORM, and the reason it is not a tuned parameter: the winning start is chosen
LEAVE-FOLD-OUT -- on the four training folds and applied to the held-out fold -- so no
deployable parameter is fitted on the native RMSD of the target it scores (contract rule 11).
The ORACLE global choice is reported beside it as the ceiling.

BASIS: built-chain CA-RMSD, tuning126.  All arms come from branch rows built in one job from
one cloud, so every contrast is paired in-job.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST                                          # noqa: E402
from s12 import instrument as I                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")


def main():
    rows = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_R_branches_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    rows[r["pdb"]] = r
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in rows]
    n = len(pdbs)
    F = np.array([rows[p]["fold"] for p in pdbs], int)
    prod = np.array([rows[p]["prod_chain"] for p in pdbs])

    tags, M = None, []
    for p in pdbs:
        fam = np.array(rows[p]["fam"])
        tg = np.array(rows[p]["tag"])
        rm = np.array(rows[p]["rmsd_nat"])
        m = fam == "GEN4D"
        o = np.argsort(tg[m])
        if tags is None:
            tags = [str(x) for x in tg[m][o]]
        assert [str(x) for x in tg[m][o]] == tags, "GEN4D columns are not aligned"
        M.append(rm[m][o])
    M = np.array(M)

    k = int(np.argmin(M.mean(0)))
    orc = ST.compare(M[:, k], prod, F, names=pdbs,
                     label="ORACLE global GEN4D start - PROD")
    lfo = np.empty(n)
    pick = {}
    for q in sorted(set(F.tolist())):
        tr = F != q
        kk = int(np.argmin(M[tr].mean(0)))
        pick[int(q)] = tags[kk]
        lfo[F == q] = M[F == q, kk]
    lc = ST.compare(lfo, prod, F, names=pdbs,
                    label="LEAVE-FOLD-OUT GEN4D start - PROD (DEPLOYABLE)")

    out = {"n": n, "INCOMPLETE": n < 126, "tags": tags,
           "per_start_mean": M.mean(0).tolist(),
           "prod_mean": float(prod.mean()),
           "best_of_4_ORACLE_mean": float(M.min(1).mean()),
           "oracle_global_start": tags[k], "oracle_global_mean": float(M[:, k].mean()),
           "oracle_vs_prod": orc, "ORACLE_NOT_DEPLOYABLE_for_oracle_arm": True,
           "lfo_pick_per_fold": pick, "lfo_mean": float(lfo.mean()), "lfo_vs_prod": lc,
           "split_half_note": "GEN4D's split-half transfer inside its own family is "
                              "-0.1338 (64% of its oracle); against PRODUCTION it is "
                              "nothing, because production's objective argmin already "
                              "reaches the same start.",
           "note": "EXPLORATORY R-16."}
    ST.save_atomic(os.path.join(RESULTS, "s32_R_gen4d_start.json"), out,
                   module_file=__file__)

    print("n = %d%s" % (n, "  (INCOMPLETE)" if n < 126 else ""))
    print("GEN4D named starts (deg): %s" % tags)
    print("per-start mean chain RMSD: %s" % np.round(M.mean(0), 4).tolist())
    print("PROD %.4f | best-of-4 ORACLE %.4f" % (prod.mean(), M.min(1).mean()))
    print("\nORACLE global start %s -> %.4f | vs PROD eff %+.4f %.2fxMDE %dW/%dL  %s"
          % (tags[k], M[:, k].mean(), orc["effect"], abs(orc["effect_over_mde"]),
             orc["n_better"], orc["n_worse"], orc["verdict"][:46]))
    print("LEAVE-FOLD-OUT pick per fold: %s" % pick)
    print("LFO arm %.4f | vs PROD eff %+.4f SE %.4f MDE %.4f %.2fxMDE %dW/%dL "
          "foldCI[%+.4f,%+.4f] %d/5"
          % (lfo.mean(), lc["effect"], lc["se"], lc["mde"], abs(lc["effect_over_mde"]),
             lc["n_better"], lc["n_worse"], lc["ci95_fold"][0], lc["ci95_fold"][1],
             lc["folds_same_sign"]))
    print("  %s" % lc["verdict"][:90])
    print("\nwrote", os.path.join(RESULTS, "s32_R_gen4d_start.json"))


if __name__ == "__main__":
    main()
