"""s22/mgeomrouter.py -- DOES THE WITHIN-WINDOW DIVERSITY PROFILE ROUTE THE m-LADDER?

The one honest attempt at the coordinator's proposed feature, after `mreal.py`'s split-half showed
per-target optimal m IS a stable, transferable property (65% of the apparent headroom survives an
independent half-pool, 2.4x its own MDE) while three native-free routers (mine, twice; the audit
lane's) captured ~0% of it using target-level distogram-confidence/score-distribution features.

METHOD, held-out-fold, exactly as `s22/routercv.py`: one Ridge regression per m-ladder arm
(`avg_500,150,75,20,5,1`), trained on the NEW `s22/mgeom.py` geometry features only (declared
primary -- no mixing with the old `rg_*`/`score_*` family, so a positive result is attributable to
the new feature class and not to smuggled old signal), alpha by generalised CV within training
folds, predicted on the held-out fold, routed by argmin.  A combined (old+new) feature set is run
alongside as a secondary, not a replacement for the declared primary.

The ORACLE best-m is used ONLY inside training (as the regression target IS each arm's own true
RMSD, which already is exactly what mreal.py calls "using half A's RMSD to native" -- there is no
separate labelling step here beyond what routercv.py already does).  The reported number is the
held-out-fold RMSD of the CHOSEN arm, never a label-accuracy score.

FALSIFIER, stated before this file's own number exists: dead unless the geometry-feature router
beats the fixed m=75 incumbent on the FULL POOL (3.048, not mreal.py's half-pool 3.072 -- the two
are not comparable, per the coordinator's caution) past its own MDE with a CI excluding zero, on
the m-ladder-only arm set.  A result that captures some but not most of mreal.py's 65% transferable
signal is reported at its own measured fraction, not rounded up to "it works" or down to "null".
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from sklearn.linear_model import RidgeCV

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s15 import seed as SD                                          # noqa: E402
from s22.routercv import ALL_FEATS, INCUMBENT, ALPHAS, _stat, _fmt  # noqa: E402
M_LADDER = ["avg_500", "avg_150", "avg_75", "avg_20", "avg_5", "avg_1"]

GEOM_FEATS = ["spread_500", "spread_150", "spread_75", "spread_20", "spread_5",
              "spread_ratio_150", "spread_ratio_75", "spread_ratio_20", "spread_ratio_5",
              "dspread_500_150", "dspread_150_75", "dspread_75_20", "dspread_20_5",
              "scorefrac_150", "scorefrac_75", "scorefrac_20", "scorefrac_5"]


def _load_joined():
    rd = json.load(open(os.path.join(RESULTS, "routerdata.json")))
    assert rd["complete"]
    mg = json.load(open(os.path.join(RESULTS, "mgeom.json")))
    assert mg["complete"]
    geo = {r["pdb"]: r for r in mg["rows"]}
    rows = []
    for r in rd["rows"]:
        g = geo[r["pdb"]]
        merged = dict(r)
        for k in GEOM_FEATS:
            merged[k] = g[k]
        rows.append(merged)
    return rows


def cv_router(rows, arm_names, feat_names, folds_arr):
    n = len(rows)
    X = np.array([[r[f] for f in feat_names] for r in rows], float)
    Y = np.array([[r["arms"][a] for a in arm_names] for r in rows], float)
    chosen_rmsd = np.full(n, np.nan)
    chosen_arm = [None] * n
    for f in sorted(set(folds_arr.tolist())):
        tr = folds_arr != f
        te = folds_arr == f
        mu, sd = X[tr].mean(0), X[tr].std(0); sd[sd == 0] = 1.0
        Xtr = (X[tr] - mu) / sd; Xte = (X[te] - mu) / sd
        preds = np.zeros((te.sum(), len(arm_names)))
        for ai in range(len(arm_names)):
            m = RidgeCV(alphas=ALPHAS)
            m.fit(Xtr, Y[tr, ai])
            preds[:, ai] = m.predict(Xte)
        pick = np.argmin(preds, axis=1)
        te_idx = np.where(te)[0]
        for k, idx in enumerate(te_idx):
            chosen_arm[idx] = arm_names[pick[k]]
            chosen_rmsd[idx] = Y[idx, pick[k]]
    assert not np.isnan(chosen_rmsd).any()
    return chosen_rmsd, chosen_arm


def main():
    rows = _load_joined()
    folds = np.array([r["fold"] for r in rows], int)
    inc = np.array([r["arms"][INCUMBENT] for r in rows], float)
    oracle_m = np.array([min(r["arms"][a] for a in M_LADDER) for r in rows], float)
    headroom = float((inc - oracle_m).mean())
    print("n=%d.  m-ladder ORACLE %.4f  incumbent(full pool) %.4f  headroom %.4f"
          % (len(rows), oracle_m.mean(), inc.mean(), headroom))
    print("(mreal.py's half-pool numbers -- fixed 3.072, selected 2.833, oracle 2.705 -- are NOT")
    print(" comparable to these full-pool numbers; not mixed into this table.)\n")

    rng = SD.stable_rng("s22", "mgeomrouter", "report")
    from collections import Counter
    for lab, feats in (("GEOMETRY ONLY (declared primary)", GEOM_FEATS),
                       ("old features (rg_*/score_*, for comparison)", ALL_FEATS),
                       ("combined (secondary)", GEOM_FEATS + ALL_FEATS)):
        chosen, arm = cv_router(rows, M_LADDER, feats, folds)
        d = chosen - inc
        s = _stat(d, rng)
        frac = 100.0 * (-s["mean"]) / headroom
        cnt = Counter(arm)
        print("%-46s achieved %.4f   %s" % (lab, chosen.mean(), _fmt(s)))
        print("   captures %.1f%% of the %.4f full-pool m-ladder headroom.  chosen: %s\n"
              % (frac, headroom, dict(cnt.most_common())))

    out = {"note": "geometry-feature m-ladder router, full-pool basis, not comparable to mreal.py's half-pool numbers"}
    path = os.path.join(RESULTS, "mgeomrouter.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(out, fh)
    os.replace(tmp, path)


if __name__ == "__main__":
    main()
