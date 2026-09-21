#!/usr/bin/env python
"""S32 LANE R -- PRODUCER for `s32/results/s32_R_sparse_control.json`.

    python s32/s32_R_sparse_analyse.py

Aggregates `s32_R_sparsectl_shard*.jsonl` (written by `s32_R_sparse_control.py`) into the
artefact S32-L(R1) quotes.

WHY IT EXISTS AS A FILE.  The aggregation was first done in an inline shell heredoc, so the
artefact had **no committed producer** -- the same defect lane V's AUDIT 11 caught on the two
dilation artefacts, found here by lane R's own producer audit before anyone else had to.
`findings-prose-is-not-evidence-of-code`: a result nobody can re-run is not a result.

BASIS: built-chain CA-RMSD, tuning126, n = 126.  `PROD` and every sparse arm are projected in
the same process per target, so every contrast is paired in-job.

WHAT THE ARTEFACT SETTLES.  The ladder's ORACLE sparse rungs project for ~free; the obvious
reading is that SPARSITY is what makes projection cheap.  `RANDSPARSE` is the matched control
-- same s = 10, same averaging operator, same projection, same job, three pinned draws --
differing from the ORACLE rung ONLY in whether the native chose the members.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST                                          # noqa: E402
from s12 import instrument as I                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
ARMS = ["PROD", "SCORESPARSE", "RANDSPARSE_d0", "RANDSPARSE_d1", "RANDSPARSE_d2"]


def main():
    R = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_R_sparsectl_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    R[r["pdb"]] = r
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in R]
    n = len(pdbs)
    F = np.array([R[p]["fold"] for p in pdbs], int)

    def col(k, fld):
        return np.array([R[p]["arms"][k][fld] for p in pdbs], float)

    out = {"n": n, "INCOMPLETE": n < 126,
           "basis": ("built chain CA-RMSD, tuning126; PROD and every arm projected in the "
                     "same process per target"),
           "prereg": ("EXPLORATORY R-12 in s32/MULTIPLICITY.md; the control S32-L9 demands. "
                      "Matched to the operator's own space: same s, same averaging operator, "
                      "same projection call, same job."),
           "arms": {}}
    print("n = %d%s\n" % (n, "  (INCOMPLETE)" if n < 126 else ""))
    print("%-16s %8s %8s %8s | %9s %9s %9s | %8s"
          % ("arm", "cloud", "chain", "d", "price", "orthog", "cos", "cos SE"))
    print("-" * 92)
    for k in ARMS:
        e, d, ch, pr, nu, co = (col(k, x) for x in
                                ("cloud_rmsd", "d_to_cloud", "chain_rmsd",
                                 "price_observed", "price_isotropic_null", "cos_align"))
        out["arms"][k] = dict(
            ORACLE_NOT_DEPLOYABLE=False,
            cloud_mean=float(e.mean()), chain_mean=float(ch.mean()),
            chain_se=float(ch.std(ddof=1) / math.sqrt(n)),
            d_mean=float(d.mean()),
            price_mean=float(pr.mean()),
            price_se=float(pr.std(ddof=1) / math.sqrt(n)),
            orthog_null_mean=float(nu.mean()),
            cos_mean=float(co.mean()), cos_se=float(co.std(ddof=1) / math.sqrt(n)),
            cos_median=float(np.median(co)), n_cos_positive=int((co > 0).sum()),
            price_minus_null=ST.compare(pr, nu, F, names=pdbs,
                                        label="%s: price - orthogonal null" % k))
        print("%-16s %8.4f %8.4f %8.4f | %+9.4f %+9.4f %+9.4f | %8.4f"
              % (k, e.mean(), ch.mean(), d.mean(), pr.mean(), nu.mean(), co.mean(),
                 co.std(ddof=1) / math.sqrt(n)))

    rs = np.array([out["arms"]["RANDSPARSE_d%d" % i]["price_mean"] for i in range(3)])
    rc = np.array([out["arms"]["RANDSPARSE_d%d" % i]["cos_mean"] for i in range(3)])
    out["RANDSPARSE_draw_distribution"] = dict(
        price_draw_means=rs.tolist(), price_draw_mean=float(rs.mean()),
        price_draw_sd=float(rs.std()), cos_draw_means=rc.tolist(),
        cos_draw_mean=float(rc.mean()), cos_draw_sd=float(rc.std()), n_draws=3,
        note="contract rule 10: a random control needs its own distribution, not its best draw")
    print("\nRANDSPARSE draw distribution (contract rule 10):")
    print("  price %s -> mean %+.4f sd %.4f" % (np.round(rs, 4).tolist(), rs.mean(), rs.std()))
    print("  cos   %s -> mean %+.4f sd %.4f" % (np.round(rc, 4).tolist(), rc.mean(), rc.std()))
    for k in ARMS:
        c = out["arms"][k]["price_minus_null"]
        print("  %-14s price-null %+.4f  %.2fxMDE  %dW/%dL  %s"
              % (k, c["effect"], abs(c["effect_over_mde"]), c["n_better"], c["n_worse"],
                 c["verdict"][:36]))
    ST.save_atomic(os.path.join(RESULTS, "s32_R_sparse_control.json"), out,
                   module_file=__file__)
    print("\nwrote", os.path.join(RESULTS, "s32_R_sparse_control.json"))


if __name__ == "__main__":
    main()
