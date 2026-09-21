#!/usr/bin/env python
"""S32 LANE R -- the CANONICAL ALIGNED VIEW of the branch rows.

WHY THIS FILE EXISTS.  `s32_R_branches_shard*.jsonl` is HETEROGENEOUS: rows written before
09:20 on 2026-09-21 carry 308 branches (three random draws) and rows after carry 158 (one).
Lane R's own analysis pins the analysed set to `GEN4 + GEN4D + MEM75 + RAND0`, but lane V
read the raw rows and got

    col   0-3    GEN4     aligned
    col   4-78   MEM75    aligned
    col  79-153  RAND0    aligned
    col 154-157  RAND1 on the 308-rows  <->  GEN4D on the 158-rows   MIS-ALIGNED

which (a) made its ORACLE-best-branch a per-target minimum over a VARIABLE K, overstating
the ceiling, and (b) biased `split_half_transfer` toward zero, because a column did not mean
the same thing on every target.  Both are lane R's fault: the shape of an artefact other
lanes were already consuming was changed mid-run.

This module emits ONE object that cannot have that problem:

    s32/results/s32_R_branches_aligned.json
        pdbs        (n,)            pinned tuning126 order, restricted to completed rows
        folds       (n,)
        prod_chain  (n,)            production's chain RMSD, from the SAME job as the branches
        columns     (158,)          "<family>|<tag>", identical for every target
        rmsd_nat    (n, 158)        ORACLE LABEL
        <criterion> (n, 158)        every native-free per-branch criterion

    Column k means the same branch construction on every target, so a column may be chosen
    on one half of the targets and scored on the other.

`typicality` IS A DISTANCE -- lower = MORE typical.  Minimise it to select the most typical
branch.  Stated here because the name already inverted one cross-lane claim.
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

from s12 import instrument as I                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
FAMILIES = ("GEN4", "MEM75", "RAND0", "GEN4D")
CRITERIA = ["rama_nlp", "rama20_nlp", "ramah", "posphi_frac", "disto_risk", "disto_mae",
            "legacy", "typicality", "obj1", "obj0", "d_to_C", "d0_to_C", "rg",
            "vbond_mean", "vbond_sd", "d_to_cloud"]


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

    cols = None
    keep_of = {}
    for p in pdbs:
        r = rows[p]
        fam = np.asarray(r["fam"]).astype(str)
        tag = np.asarray(r["tag"]).astype(str)
        idx = []
        for fm in FAMILIES:
            w = np.where(fam == fm)[0]
            idx.extend(w[np.argsort(tag[w])].tolist())      # a deterministic within-family order
        idx = np.array(idx, int)
        names = ["%s|%s" % (fam[k], tag[k]) for k in idx]
        if cols is None:
            cols = names
        elif names != cols:
            raise SystemExit("column names differ on %s -- cannot align:\n  %s\n  %s"
                             % (p, names[:4], cols[:4]))
        keep_of[p] = idx

    out = {
        "n": len(pdbs), "n_columns": len(cols), "columns": cols,
        "families": list(FAMILIES),
        "pdbs": pdbs,
        "folds": [int(rows[p]["fold"]) for p in pdbs],
        "prod_chain": [float(rows[p]["prod_chain"]) for p in pdbs],
        "cloud_rmsd": [float(rows[p]["cloud_rmsd"]) for p in pdbs],
        "basis": "built chain CA-RMSD, tuning126; branches and prod_chain from the same job",
        "ORACLE": "rmsd_nat is a LABEL; no column construction or criterion reads the native",
        "typicality_direction": "typicality is a DISTANCE to the pool: LOWER = MORE typical",
        "why": ("the raw shard rows are heterogeneous (308 branches before 09:20, 158 after) "
                "and columns 154-157 are RAND1 on the wide rows and GEN4D on the narrow "
                "ones; this view is aligned by (family, tag) so a column means the same "
                "construction on every target"),
        "rmsd_nat": [np.asarray(rows[p]["rmsd_nat"], float)[keep_of[p]].tolist()
                     for p in pdbs],
    }
    for c in CRITERIA:
        if c not in rows[pdbs[0]]:
            continue
        out[c] = [np.asarray(rows[p][c], float)[keep_of[p]].tolist() for p in pdbs]

    path = os.path.join(RESULTS, "s32_R_branches_aligned.json")
    tmp = path + ".%d.tmp" % os.getpid()            # a SHARED .tmp is not atomic
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(out, fh)
    os.replace(tmp, path)
    print("n = %d targets x %d columns" % (out["n"], out["n_columns"]))
    print("families: %s" % (FAMILIES,))
    print("first columns: %s" % cols[:3])
    print("cols 154-157:  %s" % cols[154:158])
    print("criteria present: %s" % [c for c in CRITERIA if c in out])
    print("wrote", path)


if __name__ == "__main__":
    main()
