#!/usr/bin/env python
"""s26/w_provenance_readout.py -- H_P3 of s26/PREREG_window_provenance.md: the class-weighted
readout.  The shipped top-75 (production `sub`) averaged in the production medoid frame with
per-member weights by class (fragment-class relative weight w in {0, 0.5, 1, 2}; peptide-derived
members weight 1), one projection per weight; the weight chosen leave-fold-out on the BUILT
CHAIN; the control: the same weights permuted across the 75 members (seeded draws) at the
chosen weight.  GATED (the LFO choice reads natives); the clouds and chains of the four weight
arms are computed first and stored, then the natives are read.

Declared reduction (stated as a deviation in the ledger): the permuted control uses 4 draws per
target at the chosen weight instead of the prereg's 8, to keep the job under one hour on a
loaded box (126 x (4 + 4) projections).

    python s26/jobrun.py --agent W --tag CPU --name w_provenance_readout --est-ram 0.3 -- python s26/w_provenance_readout.py
"""
from __future__ import annotations

import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
import w_selfcopy as W                               # noqa: E402

K, TOPM = W.K, W.TOPM
WEIGHTS = (0.0, 0.5, 1.0, 2.0)
N_PERM = 4
OUT = os.path.join(W.RES, "w_selfcopy_provenance_readout.json")


def weighted_cloud(top, w, ref):
    """Weighted mean of the members superposed onto the production medoid."""
    M = I.superpose_batch(np.asarray(top, float), np.asarray(ref, float))
    w = np.asarray(w, float)
    if w.sum() <= 0:
        return None
    return (M * w[:, None, None]).sum(0) / w.sum()


def main():
    W.require_signoff("provenance readout")
    from s15 import seed as SD
    Z = W.load_result("provenance_census")
    assert Z and Z.get("complete")
    tg = W.targets_by_pdb()
    prev = json.load(open(OUT)) if os.path.exists(OUT) else None
    rows = prev["rows"] if prev else []
    done = {r["pdb"] for r in rows}
    t_last = time.time()
    for r in Z["rows"]:
        pdb = r["pdb"]
        if pdb in done:
            continue
        t = tg[pdb]; u = I.load_univ(pdb); N = np.asarray(u["nat_ca"], float)
        order = np.asarray(u["order"], int); pool = order[:K]
        rec = I.shipped_record(pdb); sub = np.asarray(rec["sub"], int)
        top = np.asarray(u["W"][pool][sub], float)
        cls = np.array(r["class_of_pool_member"])[sub]
        is_frag = (cls == "fragment")
        P = I.pairwise_rmsd(top); ref = top[I.medoid(P)]
        row = {"pdb": pdb, "fold": t["fold"], "n_frag_in_top75": int(is_frag.sum()), "arms": {}, "perm": {}}
        for wf in WEIGHTS:
            w = np.where(is_frag, wf, 1.0)
            C = weighted_cloud(top, w, ref)
            if C is None:                                # no peptide member and wf = 0: undefined, fall back to uniform
                C = weighted_cloud(top, np.ones(TOPM), ref); fallback = True
            else:
                fallback = False
            pr = I.project(C, t["seq"], int(t["fold"]))
            row["arms"]["%g" % wf] = {"cloud": float(I.ca_rmsd(C, N)), "arm": float(I.ca_rmsd(pr["ca"], N)), "fit": float(I.ca_rmsd(pr["fit_ca"], N)), "fallback": fallback}
            if wf == 1.0:
                assert abs(np.abs(C - np.asarray(rec["avg_ca"], float)).max()) < 1e-9, "uniform weights do not reproduce avg_ca"
        # permuted-weight control at every non-uniform weight (chosen weight is picked later, LFO)
        for wf in WEIGHTS:
            if wf == 1.0:
                continue
            for q in range(N_PERM):
                rng = SD.stable_rng("w_provenance_readout", pdb, wf, q)
                w = np.where(is_frag, wf, 1.0)[rng.permutation(TOPM)]
                C = weighted_cloud(top, w, ref)
                if C is None:
                    C = weighted_cloud(top, np.ones(TOPM), ref)
                pr = I.project(C, t["seq"], int(t["fold"]))
                row["perm"]["%g:%d" % (wf, q)] = {"cloud": float(I.ca_rmsd(C, N)), "arm": float(I.ca_rmsd(pr["ca"], N))}
        rows.append(row)
        print("  [%3d/126] %s frag %2d/75: arm w0 %.3f w0.5 %.3f w1 %.3f w2 %.3f" % (len(rows), pdb, row["n_frag_in_top75"], row["arms"]["0"]["arm"], row["arms"]["0.5"]["arm"], row["arms"]["1"]["arm"], row["arms"]["2"]["arm"]), flush=True)
        if time.time() - t_last > 300 or len(rows) % 10 == 0:
            W.save("provenance_readout", {"label": "H_P3 class-weighted readout, GATED", "rows": rows}, rows=rows, n_expected=len(Z["rows"]))
            t_last = time.time()
        del u
    pdbs = [r["pdb"] for r in rows]; folds = ST.pinned_folds(pdbs)
    A = {wf: np.array([r["arms"]["%g" % wf]["arm"] for r in rows]) for wf in WEIGHTS}
    Cc = {wf: np.array([r["arms"]["%g" % wf]["cloud"] for r in rows]) for wf in WEIGHTS}
    # leave-fold-out choice of the fragment weight on the built chain
    choice = np.empty(len(rows), dtype=float)
    for f in sorted(set(folds.tolist())):
        train = folds != f
        means = [A[wf][train].mean() for wf in WEIGHTS]
        choice[folds == f] = WEIGHTS[int(np.argmin(means))]       # ties: first in grid order
    routed_arm = np.array([r["arms"]["%g" % choice[k]]["arm"] for k, r in enumerate(rows)])
    routed_cloud = np.array([r["arms"]["%g" % choice[k]]["cloud"] for k, r in enumerate(rows)])
    perm_arm = np.array([np.mean([r["perm"]["%g:%d" % (choice[k], q)]["arm"] for q in range(N_PERM)]) if choice[k] != 1.0 else r["arms"]["1"]["arm"] for k, r in enumerate(rows)])
    perm_cloud = np.array([np.mean([r["perm"]["%g:%d" % (choice[k], q)]["cloud"] for q in range(N_PERM)]) if choice[k] != 1.0 else r["arms"]["1"]["cloud"] for k, r in enumerate(rows)])
    out = {"label": "H_P3 class-weighted readout, GATED", "rows": rows, "prereg": "s26/PREREG_window_provenance.md",
           "chosen_weight_by_fold": {int(f): float(choice[folds == f][0]) for f in sorted(set(folds.tolist()))}, "stats": {}}
    for name, a, b, basis in (("routed_minus_uniform", routed_arm, A[1.0], "arm"), ("routed_minus_permuted", routed_arm, perm_arm, "arm"),
                              ("routed_minus_uniform", routed_cloud, Cc[1.0], "cloud"), ("routed_minus_permuted", routed_cloud, perm_cloud, "cloud")):
        res = ST.compare(a, b, folds, names=pdbs, label="%s [%s], fragment weight chosen leave-fold-out on the built chain" % (name, basis))
        out["stats"]["%s_%s" % (name, basis)] = res; print(ST.fmt(res))
    for wf in WEIGHTS:
        if wf == 1.0:
            continue
        res = ST.compare(A[wf], A[1.0], folds, names=pdbs, label="fixed fragment weight %g minus uniform [arm] (all folds; NOT the routed arm)" % wf)
        out["stats"]["fixed_w%g_minus_uniform_arm" % wf] = res; print(ST.fmt(res))
    print("  chosen weight by fold:", out["chosen_weight_by_fold"])
    W.save("provenance_readout", out, rows=rows, n_expected=len(Z["rows"]), complete_keys=("pdb", "arms", "perm"))
    print("  ->", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
