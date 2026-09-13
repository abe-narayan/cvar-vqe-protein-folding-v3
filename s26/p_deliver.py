"""s26/p_deliver.py -- hand the best C2 rung's BUILT CHAINS to lane PH for C3 stage 2.

Format (coordinator, 2026-09-13 09:00): `s26/results/p_best_rung_chains.json`, written with
`ST.save_atomic(..., complete_keys=..., rows=rows, n_expected=126, module_file=__file__)`; rows
keyed by pdb (a `rows` list for the completeness gate plus an `index` pdb -> row position), each
row carrying `pdb`, `rung`, `phi` and `psi` (the built chain's torsions in RADIANS exactly as
`s12.instrument.project` -> `core.project.lam_path` emits them at lam = 0.3, lists of length n),
`ca` (the built-chain CA coordinates, n x 3) and `rmsd_arm` (the rung's built-chain CA-RMSD).

Why a separate module: `s26/p_ladder.py` is the source of the training jobs in flight and the
contract forbids editing a module while a job launched from it is running.  This module imports
it unchanged and re-runs the rung's deployable path plus one projection per target (2.9 s), then
asserts that the re-projected chain's RMSD equals the value in the rung's eval JSON, so the file
is the same emission the verdict was read from.  Gated on "PHASE 0 SIGNED OFF".

    python s26/p_deliver.py --rung <rung> [--seed 0]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from s26 import p_ladder as L                        # noqa: E402

OUT = os.path.join(HERE, "results", "p_best_rung_chains.json")
NEED = ["pdb", "rung", "n", "fold", "seq", "phi", "psi", "ca", "rmsd_arm", "rmsd_cloud", "rmsd_fit"]


def deliver(rung, seed=0, lam=None, verbose=True):
    if not L._signed_off():
        raise SystemExit("PHASE GATE: delivery reads native RMSDs; refused before sign-off.")
    L.install_guard()
    ev = json.load(open(L.result_path(rung, seed)))
    assert ev.get("complete"), "eval of rung %s is not complete" % rung
    evrows = {r["pdb"]: r for r in ev["rows"]}
    models = L.train_rung(rung, seed=seed, verbose=False) if rung in L.TRAINABLE else None
    rows = []
    for c, t in enumerate(I.targets()):
        pdb = t["pdb"]; u = I.load_univ(pdb)
        idx = I.pool_idx(u); W_pool = np.asarray(u["W"][idx], float)
        if rung == "mix":
            P0, i, j = L.posterior("shipped", pdb, t["seq"], t["fold"])
            Dp = I.pair_dists(W_pool, i, j)
            Q, i, j = L.posterior("mix", pdb, t["seq"], t["fold"], lam=lam, D_pool=Dp)
        else:
            Q, i, j = L.posterior(rung, pdb, t["seq"], t["fold"], models=models, seed=seed)
        sc, o, C, rt, Dp = L.score_target(Q, i, j, t["seq"], W_pool)
        pr = I.project(C, t["seq"], int(t["fold"]))
        nat = np.asarray(u["nat_ca"], float)
        arm = float(I.ca_rmsd(pr["ca"], nat))
        key = "arm" if rung != "mix" else "arm@%g" % lam
        ref = float(evrows[pdb][key])
        assert abs(arm - ref) < 1e-6, "re-projection differs from the eval JSON on %s: %.6f vs %.6f" % (pdb, arm, ref)
        rows.append({"pdb": pdb, "rung": rung if rung != "mix" else "mix@%g" % lam, "n": int(t["n"]), "fold": int(t["fold"]),
                     "seq": t["seq"], "phi": [float(x) for x in pr["phi"]], "psi": [float(x) for x in pr["psi"]],
                     "ca": [[float(v) for v in row] for row in pr["ca"]], "rmsd_arm": arm,
                     "rmsd_cloud": float(I.ca_rmsd(C, nat)), "rmsd_fit": float(I.ca_rmsd(pr["fit_ca"], nat))})
        del u
        if verbose and (c + 1) % 25 == 0:
            print("  %d/126" % (c + 1), flush=True)
    out = {"rung": rung, "seed": seed, "lam": lam, "torsion_units": "radians (core.project.lam_path emission, lam=0.3)",
           "rows": rows, "index": {r["pdb"]: k for k, r in enumerate(rows)},
           "mean_rmsd_arm": float(np.mean([r["rmsd_arm"] for r in rows])),
           "eval_artefact": L.result_path(rung, seed)}
    ST.save_atomic(OUT, out, complete_keys=NEED, rows=rows, n_expected=126, module_file=__file__)
    print("wrote %s: rung %s, mean built-chain RMSD %.4f" % (OUT, out["rung"], out["mean_rmsd_arm"]))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rung", required=True); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--lam", type=float, default=None, help="mix only: the leave-fold-out lam to deliver")
    a = ap.parse_args()
    deliver(a.rung, a.seed, a.lam)
