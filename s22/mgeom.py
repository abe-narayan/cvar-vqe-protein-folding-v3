"""s22/mgeom.py -- THE WITHIN-WINDOW DIVERSITY PROFILE, the coordinator's proposed feature for the
m-ladder router, after mreal.py's split-half REFUTED the "m*-headroom is pure min-of-K noise"
hypothesis: 65% of the apparent headroom TRANSFERS across independent halves of the same target's
pool (-0.239 [-0.309,-0.173], 2.4x its own MDE) -- the optimum m is a real, stable per-target
property.  Three native-free routers (my Ridge on 5 arm sets, my RF in-fold/held-out, the audit
lane's length router) have all captured ~0% of it.  This file tests ONE new, mechanistically
distinct feature family before concluding the blindness is structural rather than a feature gap.

WHY THIS FEATURE, not another target-level summary.  What `m` controls is how far down the
score-ordering the AVERAGING BENEFIT persists -- i.e. how the retained set's own internal geometry
changes as the window widens.  That is a property of the CANDIDATE SET, not of the distogram's
confidence in its own prediction (which was already tried in `s22/routerdata.py`'s `score_*` family
and in the coordinator's own distogram-confidence test, both null on which-m).  The natural
observable is:

    spread(m) = mean pairwise RMSD among the top-m candidates BY SCORE, m in (500,150,75,20,5)

and its SHAPE: ratios to the full-pool spread, and the successive differences (where it knees).
Also carried: the SCORE profile itself at the same rungs, `scorefrac(m)`, since a set whose score
gap to the boundary is large is a qualitatively different selection event from one whose scores are
smooth there -- both are candidate-set-internal, native-free, and computed from data the pipeline
already retrieves.

Declared BEFORE any m-router number exists on these features: `m=1` carries no internal pairwise
spread by construction and is excluded from the spread family; `spread(500)` is the full-pool
value and the natural normaliser.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402

MS_SPREAD = (500, 150, 75, 20, 5)


def _save(obj):
    path = os.path.join(RESULTS, "mgeom.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def run():
    tg = I.targets()
    rows = []
    t0 = time.time()
    print("targets: %d" % len(tg), flush=True)
    for c, t in enumerate(tg):
        pdb, n = t["pdb"], int(t["n"])
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(n)
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        order = np.argsort(sc, kind="stable")
        Ws = W[order]
        scs = sc[order]
        P = I.pairwise_rmsd(Ws)          # (500,500), sorted by score -- top-m is P[:m,:m]

        row = {"pdb": pdb, "n": n, "fold": int(t["fold"])}
        spreads = {}
        for m in MS_SPREAD:
            sub = P[:m, :m]
            spreads[m] = float(sub[np.triu_indices(m, 1)].mean())
            row["spread_%d" % m] = spreads[m]
        s500 = spreads[500]
        for m in (150, 75, 20, 5):
            row["spread_ratio_%d" % m] = spreads[m] / s500 if s500 > 0 else 0.0
        seq_ms = (500, 150, 75, 20, 5)
        for a, b in zip(seq_ms[:-1], seq_ms[1:]):
            row["dspread_%d_%d" % (a, b)] = spreads[a] - spreads[b]
        rng_sc = scs[499] - scs[0]
        for m in (150, 75, 20, 5):
            row["scorefrac_%d" % m] = float((scs[m - 1] - scs[0]) / rng_sc) if rng_sc > 0 else 0.0
        rows.append(row)
        if (c + 1) % 10 == 0:
            print("  %d/%d  (%.1fs)" % (c + 1, len(tg), time.time() - t0), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    ok = len(rows) == len(tg)
    feat_keys = sorted(k for k in rows[0] if k not in ("pdb", "n", "fold"))
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "feat_keys": feat_keys})
    print("DONE complete=%s" % ok)
    print("features:", feat_keys)
    return rows


if __name__ == "__main__":
    run()
