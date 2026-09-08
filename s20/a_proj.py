"""SPRINT 20, AGENT A -- the BUILT-CHAIN pass.  Phase 1 emits point clouds; this emits structures.

`I.coordinate_average` returns a contracted point cloud (mean virtual Ca-Ca bond 2.961 A against
the physical 3.804 A).  Every number in `a_src.py` is therefore in the POINT CLOUD basis and is
comparable only to the 3.048 A point-cloud reference.  This module puts every arm on the
ideal-geometry manifold with the pipeline's own operator --

    X_built = I.project(X_cloud, seq, fold)["fit_ca"]        the lam = 0 arm

-- which is the incumbent's own basis, mean 3.204 A.  A source that only looked good by
contracting the chain loses that here, by construction.

Reads the candidate sets `a_src.py` cached in `s20/cache/sets_<pdb>.npz`, so nothing is rescored
and the two phases consume bit-identical candidate sets.

Run:  python -m s20.a_proj            (checkpoints every target; resumable)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I              # noqa: E402
from s20 import a_src as S                   # noqa: E402
from s20.a_gate import _partial              # noqa: E402

OUT = os.path.join(S.RESULTS, "a_proj.json")

SRC = ["pool", "pep", "prot", "halfA", "halfB", "rand500", "tors", "unsel"]
FUSE = [("pool", "pep"), ("halfA", "halfB"), ("pool", "tors"), ("pool", "unsel")]
MERGE = [("pool", "pep"), ("halfA", "halfB")]
# rho on the BUILT basis, which is the basis the falsifier should be read on: every built chain
# has the same ideal virtual bond (3.804 A), so no arm's alignment can be an artefact of a
# differently contracted point cloud.  `helix` is the zero-information reference to partial out.
RHO_REF = "pool"


def _bond(x):
    b = np.linalg.norm(np.asarray(x)[1:] - np.asarray(x)[:-1], axis=1)
    return float(b.mean()), float(b.min())


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    rows = json.load(open(OUT))["rows"] if os.path.exists(OUT) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        if pdb in done:
            continue
        z = np.load(os.path.join(S.CACHE, f"sets_{pdb}.npz"))
        u = I.load_univ(pdb)
        nat = u["nat_ca"]
        W_all = np.asarray(u["W"], float)
        e = {"pdb": pdb, "n": n, "fold": fold}

        WT, X = {}, {}
        for s in SRC:
            if s == "tors":
                WT[s] = np.asarray(I.build_ca(z["tors_phi"], z["tors_psi"]), float)
            else:
                WT[s] = W_all[np.asarray(z[f"top_{s}"], int)]
            X[s], _b = I.coordinate_average(WT[s])
        for a, b in FUSE:
            xb = I.superpose_batch(X[b][None], X[a])[0]
            X[f"fuse|{a}|{b}"] = 0.5 * (X[a] + xb)
        for a, b in MERGE:
            X[f"merge|{a}|{b}"], _ = I.coordinate_average(np.concatenate([WT[a], WT[b]], 0))

        B = {}
        for k, xc in X.items():
            pr = I.project(np.asarray(xc, float), seq, fold)
            bm, bn = _bond(pr["fit_ca"])
            e[f"rmsdP|{k}"] = float(I.ca_rmsd(pr["fit_ca"], nat))
            e[f"rmsdL|{k}"] = float(I.ca_rmsd(pr["ca"], nat))
            e[f"bondP|{k}"] = bm
            e[f"bondPmin|{k}"] = bn
            e[f"rmsdC|{k}"] = float(I.ca_rmsd(xc, nat))        # the cloud it came from
            B[k] = np.asarray(pr["fit_ca"], float)

        # ---- the falsifier metric, on the BUILT basis
        i, j = I.pair_index(n)
        dtrue = np.linalg.norm(nat[i] - nat[j], axis=1)
        hel = I.build_ca(np.full(n, np.deg2rad(-63.0)), np.full(n, np.deg2rad(-42.0)))
        e_hel = np.linalg.norm(hel[i] - hel[j], axis=1) - dtrue
        EB = {k: np.linalg.norm(v[i] - v[j], axis=1) - dtrue for k, v in B.items()}
        EB["helix"] = e_hel
        nm = list(EB)
        Z = S._std(np.stack([EB[k] for k in nm]))
        C = Z @ Z.T
        r0 = nm.index(RHO_REF)
        for a_ in range(len(nm)):
            e[f"rhoB|{nm[a_]}"] = float(C[r0, a_])
            if nm[a_] != "helix":
                e[f"rhoBP|{nm[a_]}"] = _partial(EB[RHO_REF], EB[nm[a_]], e_hel)
                e[f"xrmsdB|{nm[a_]}"] = float(I.ca_rmsd(B[RHO_REF], B[nm[a_]]))
        for a_, b_ in [("pep", "prot"), ("halfA", "halfB")]:
            ia, ib = nm.index(a_), nm.index(b_)
            e[f"rhoB|{a_}|{b_}"] = float(C[ia, ib])
            e[f"rhoBP|{a_}|{b_}"] = _partial(EB[a_], EB[b_], e_hel)
        rows.append(e)
        json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
        el = time.time() - t0
        print(f"  {len(rows)}/{len(tg)} {pdb} ({el:.0f}s, {el/max(len(rows)-len(done),1):.1f}s/tgt)",
              flush=True)
    if len(rows) == len(tg):
        open(os.path.join(S.RESULTS, "a_proj.COMPLETE"), "w").write("ok\n")
    return rows


if __name__ == "__main__":
    run()
