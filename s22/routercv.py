"""s22/routercv.py -- THE ROUTER, SCORED ONLY ON HELD-OUT FOLDS (PREREG_C.md Experiment 1).

Reads s22/results/routerdata.json (features + every candidate arm's RMSD, n=126, pinned 5-fold
structure).  Fits ONE Ridge regression per arm, per held-out fold, on the OTHER four folds only;
predicts every arm's RMSD on the held-out fold; routes each held-out target to the arm with the
lowest PREDICTED RMSD.  The achieved router number is the concatenation of five held-out
predictions -- never a number computed on the folds that trained it.

Reports, for each of five nested candidate-arm sets (so the coordinator's own family decomposition
-- s22/LEDGER.md's message on m-ladder/pool/latent -- can be cross-checked against an independently
built table):

    A          m-ladder only (6 arms, one source, one operator)
    A+medoid   all POOL arms (m-ladder + 2 consensus medoids)                    == "World A", 8 arms
    A+latent   World A + 4 latent arms (World B)                                 12 arms
    A+blend    World A + the NEW blend ladder (4 new points; blend_100==avg_75)  12 arms
    ALL        World A + World B + blend                                        16 arms

and, for each: the ORACLE ceiling (min over arms, never achievable), the best-FIXED-arm baseline
(always avg_75 = 3.048, the correct baseline per BRIEF), the matched-count RANDOM-ARM control
(sanity: routing must beat picking blind), and the ACHIEVED held-out-fold router.

Also runs the rg-only (2-feature) vs full-feature ablation on the ALL set, and reports the
descriptive per-target chosen-arm distribution by headroom tercile (Experiment 3, descriptive only,
no new falsifiable claim).
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

sys.path.insert(0, ROOT)
from s15 import seed as SD                  # noqa: E402

INCUMBENT = "avg_75"
M_LADDER = ["avg_500", "avg_150", "avg_75", "avg_20", "avg_5", "avg_1"]
MEDOIDS = ["medoid_all", "medoid_75"]
LATENT = ["lat_argmin", "lat_avg75", "lat_medoid", "lat_avg75_rand"]
BLEND_NEW = ["blend_00", "blend_25", "blend_50", "blend_75"]     # blend_100 == avg_75, excluded as dup

SETS = {
    "m_ladder_only": M_LADDER,
    "pool_all (World A)": M_LADDER + MEDOIDS,
    "A_plus_latent": M_LADDER + MEDOIDS + LATENT,
    "A_plus_blend": M_LADDER + MEDOIDS + BLEND_NEW,
    "ALL": M_LADDER + MEDOIDS + LATENT + BLEND_NEW,
}

ALL_FEATS = ["n", "rg_disto", "rg_pool_mean", "rg_pool_sd", "rg_gap", "rg_z",
             "score_mean", "score_sd", "score_gap01", "score_gap_boundary",
             "score_iqr_over_range", "score_skew", "pool_spread",
             "sim_mean_top75", "sim_mean_pool"]
RG_FEATS = ["rg_z", "rg_gap"]

ALPHAS = np.logspace(-2, 4, 40)


def _load():
    d = json.load(open(os.path.join(RESULTS, "routerdata.json")))
    assert d["complete"], "routerdata.json is not complete"
    rows = d["rows"]
    return rows


def _stat(d, rng, B=5000):
    d = np.asarray(d, float); k = len(d)
    se = float(d.std(ddof=1) / np.sqrt(k))
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return dict(mean=float(d.mean()), se=se, mde=2.8016 * se,
                lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)),
                w=int((d < 0).sum()), l=int((d > 0).sum()), n=k)


def _fmt(s):
    return "%+.4f SE %.4f MDE %.4f [%+.4f,%+.4f] %dW/%dL" % (s["mean"], s["se"], s["mde"], s["lo"], s["hi"], s["w"], s["l"])


def cv_router(rows, arm_names, feat_names, folds_arr):
    """Held-out-fold routing.  Returns per-target chosen arm and its actual RMSD, out-of-fold only."""
    n = len(rows)
    X = np.array([[r[f] for f in feat_names] for r in rows], float)
    Y = np.array([[r["arms"][a] for a in arm_names] for r in rows], float)   # (n, n_arms)
    chosen_rmsd = np.full(n, np.nan)
    chosen_arm = [None] * n
    for f in sorted(set(folds_arr.tolist())):
        tr = folds_arr != f
        te = folds_arr == f
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1.0
        Xtr = (X[tr] - mu) / sd
        Xte = (X[te] - mu) / sd
        preds = np.zeros((te.sum(), len(arm_names)))
        for ai in range(len(arm_names)):
            model = RidgeCV(alphas=ALPHAS)
            model.fit(Xtr, Y[tr, ai])
            preds[:, ai] = model.predict(Xte)
        pick = np.argmin(preds, axis=1)
        te_idx = np.where(te)[0]
        for k, idx in enumerate(te_idx):
            chosen_arm[idx] = arm_names[pick[k]]
            chosen_rmsd[idx] = Y[idx, pick[k]]
    assert not np.isnan(chosen_rmsd).any()
    return chosen_rmsd, chosen_arm


def random_control(rows, arm_names, rng, reps=500):
    n = len(rows)
    Y = np.array([[r["arms"][a] for a in arm_names] for r in rows], float)
    picks = rng.integers(0, len(arm_names), size=(reps, n))
    vals = Y[np.arange(n)[None, :], picks]     # (reps, n)
    return vals.mean(1)          # per-rep mean RMSD; report distribution


def oracle_and_fixed(rows, arm_names):
    Y = np.array([[r["arms"][a] for a in arm_names] for r in rows], float)
    inc = np.array([r["arms"][INCUMBENT] for r in rows], float)
    oracle = Y.min(1)
    winner_counts = {a: int((Y.argmin(1) == i).sum()) for i, a in enumerate(arm_names)}
    return oracle, inc, winner_counts


def main():
    rows = _load()
    n = len(rows)
    folds = np.array([r["fold"] for r in rows], int)
    print("n = %d targets, %d folds." % (n, len(set(folds.tolist()))))
    inc = np.array([r["arms"][INCUMBENT] for r in rows], float)
    print("\nBaseline (best FIXED arm, the incumbent avg_75): %.4f\n" % inc.mean())

    rng = SD.stable_rng("s22", "routercv", "report")
    results = {}
    for label, arms in SETS.items():
        oracle, inc2, wc = oracle_and_fixed(rows, arms)
        assert np.allclose(inc2, inc)
        headroom = float((inc - oracle).mean())
        rand_vals = random_control(rows, arms, rng)
        print("=" * 100)
        print("ARM SET: %-22s (%d arms): %s" % (label, len(arms), arms))
        print("  ORACLE ceiling %.4f   incumbent %.4f   HEADROOM %.4f" % (oracle.mean(), inc.mean(), headroom))
        print("  winner counts: %s" % {k: v for k, v in sorted(wc.items(), key=lambda kv: -kv[1])})
        print("  random-arm control (matched count, %d arms, 500 reps): mean %.4f  sd-over-reps %.4f"
              % (len(arms), rand_vals.mean(), rand_vals.std()))

        chosen_rmsd, chosen_arm = cv_router(rows, arms, ALL_FEATS, folds)
        d = chosen_rmsd - inc
        s = _stat(d, rng)
        per_fold = {int(f): float(d[folds == f].mean()) for f in sorted(set(folds.tolist()))}
        frac = 100.0 * (-s["mean"]) / headroom if headroom > 0 else float("nan")
        print("  ACHIEVED CV ROUTER (full features): %.4f   %s" % (chosen_rmsd.mean(), _fmt(s)))
        print("    captures %.1f%% of the %.4f headroom.  per-fold: %s" % (frac, headroom, per_fold))
        from collections import Counter
        cnt = Counter(chosen_arm)
        print("    router's chosen-arm distribution: %s" % dict(cnt.most_common()))

        results[label] = dict(arms=arms, oracle=float(oracle.mean()), incumbent=float(inc.mean()),
                               headroom=headroom, winner_counts=wc,
                               random_control_mean=float(rand_vals.mean()),
                               random_control_sd=float(rand_vals.std()),
                               achieved_router_mean=float(chosen_rmsd.mean()),
                               achieved_vs_incumbent=s, per_fold=per_fold,
                               capture_frac_pct=frac, chosen_arm_counts=dict(cnt))

    # ---------------------------------------------------------------- rg-only vs full-feature, on ALL
    print("\n" + "=" * 100)
    print("ABLATION on ALL (%d arms): rg-only (2 feats) vs full (%d feats)" % (len(SETS["ALL"]), len(ALL_FEATS)))
    for lab, feats in (("rg-only", RG_FEATS), ("full", ALL_FEATS)):
        chosen_rmsd, chosen_arm = cv_router(rows, SETS["ALL"], feats, folds)
        d = chosen_rmsd - inc
        s = _stat(d, rng)
        print("  %-10s %.4f   %s" % (lab, chosen_rmsd.mean(), _fmt(s)))
        results.setdefault("ablation", {})[lab] = dict(mean=float(chosen_rmsd.mean()), stat=s)

    # ---------------------------------------------------------------- descriptive: chosen arm by headroom tercile
    oracle_all, _, _ = oracle_and_fixed(rows, SETS["ALL"])
    headroom_pt = inc - oracle_all
    terc = np.quantile(headroom_pt, [1 / 3, 2 / 3])
    tier = np.digitize(headroom_pt, terc)
    chosen_rmsd_all, chosen_arm_all = cv_router(rows, SETS["ALL"], ALL_FEATS, folds)
    print("\n" + "=" * 100)
    print("DESCRIPTIVE (Experiment 3): CV-router's chosen arm family by per-target headroom tercile")
    from collections import Counter
    fam_of = {}
    for a in M_LADDER:
        fam_of[a] = "m_ladder"
    for a in MEDOIDS:
        fam_of[a] = "medoid"
    for a in LATENT:
        fam_of[a] = "latent"
    for a in BLEND_NEW:
        fam_of[a] = "blend"
    tercile_fams = {}
    for tlab, tval in (("low headroom", 0), ("mid headroom", 1), ("high headroom", 2)):
        idx = tier == tval
        fams = Counter(fam_of[chosen_arm_all[i]] for i in np.where(idx)[0])
        print("  %-14s (n=%3d, mean headroom %.3f): %s" % (tlab, idx.sum(), headroom_pt[idx].mean(), dict(fams)))
        tercile_fams[tlab] = dict(fams)

    out = {"results": results, "tercile_fams": tercile_fams,
           "n_targets": n, "feat_names": ALL_FEATS, "complete": True}
    path = os.path.join(RESULTS, "routercv.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)
    print("\nwrote", path)


if __name__ == "__main__":
    main()
