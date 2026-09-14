"""s26/p_rung_report.py -- the ledger block for ONE rung of the C2 ladder against the shipped posterior.

Prints ST.fmt verbatim for the built chain (PRIMARY), the point cloud and the selection endpoint,
the gamma-equivalent (gam_eff and cos in probability space and in location space, with the S25 L12
caveat), and, for the shipped rung, the anchor check against the production cache.  For `mix` the
leave-fold-out lam is chosen here (on `arm`, within the training folds) and the per-target oracle
over lam is scored against ST.best_of_k_within.  Writes s26/results/p_ladder_report_<rung>_s<seed>.json.

    python s26/p_rung_report.py --rung noesm [--seed 0]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from s26 import p_ladder as L                        # noqa: E402

PINNED = {"sel": 3.4540, "cloud": 3.0483, "arm": 3.2148, "fit": 3.2041}
CACHE_KEY = {"sel": "shipped", "cloud": "rmsd_avg", "arm": "rmsd_arm", "fit": "rmsd_fit"}


def mix_lfo(rows, folds, lams, key="arm"):
    """Leave-fold-out lam: for each fold, the lam minimising the mean `key` over the OTHER folds."""
    A = np.array([[r["%s@%g" % (key, l)] for l in lams] for r in rows])
    chosen = {}; out = np.zeros(len(rows))
    for f in sorted(set(folds.tolist())):
        tr = folds != f
        k = int(np.argmin(A[tr].mean(0)))
        chosen[int(f)] = float(lams[k]); out[folds == f] = A[folds == f, k]
    return out, chosen, A


def report(rung, seed=0):
    base = L._rows("shipped", seed)
    if base is None:
        raise SystemExit("eval --rung shipped has not completed")
    pdbs = [r["pdb"] for r in base]; folds = ST.pinned_folds(pdbs)
    out = {"rung": rung, "seed": seed, "n": len(pdbs)}
    lines = []
    if rung == "shipped":
        recs = {}
        for f in glob.glob(os.path.join(ROOT, "bench_results", "cache", I.PROD_KEY, "*.json")):
            z = json.load(open(f)); recs[z["pdb"]] = z
        lines.append("SHIPPED THROUGH THE LADDER'S OWN PATH (the anchor), n=%d" % len(pdbs))
        for k in ("sel", "cloud", "arm", "fit"):
            v = np.array([r[k] for r in base]); c = np.array([recs[p][CACHE_KEY[k]] for p in pdbs])
            d = np.abs(v - c)
            lines.append("  %-5s mean %.4f (pinned %.4f)  vs production cache per target: max abs %.2e, mean abs %.2e, n(diff>1e-6)=%d"
                         % (k, v.mean(), PINNED[k], d.max(), d.mean(), int((d > 1e-6).sum())))
            out[k] = {"mean": float(v.mean()), "pinned": PINNED[k], "maxabs_vs_cache": float(d.max()), "n_diff": int((d > 1e-6).sum())}
        g = {q: float(np.nanmean([r[q] for r in base])) for q in ("gam_prob", "cos_prob", "gam_loc", "cos_loc", "mae")}
        lines.append("  gam_eff/cos of the shipped posterior against itself: prob %+.4f/%.3f loc %+.4f/%.3f (must be 0/nan)  MAE %.4f" % (g["gam_prob"], g["cos_prob"], g["gam_loc"], g["cos_loc"], g["mae"]))
        out["progress"] = g
    else:
        rows = L._rows(rung, seed)
        if rows is None:
            raise SystemExit("eval --rung %s has not completed" % rung)
        m = {r["pdb"]: r for r in rows}; rr = [m[p] for p in pdbs]
        vals = {}
        if rung == "mix":
            lams = json.load(open(L.result_path(rung, seed)))["lams"]
            for k in ("arm", "cloud", "sel"):
                vals[k], chosen, A = mix_lfo(rr, folds, lams, k)
                if k == "arm":
                    out["lam_chosen_per_fold"] = chosen
                    bok = ST.best_of_k_within(A, n_boot=300)
                    out["oracle_over_lam"] = {q: v for q, v in bok.items() if q != "argmin_counts"}
                    lines.append("MIX: leave-fold-out lam per fold %s ; per-target ORACLE over %d lams: gain %+.4f, valid null %+.4f (share %.2f), split-half %+.4f (%.0f%% of oracle), k_eff %.2f"
                                 % (chosen, len(lams), bok["observed_gain"], bok["null_across_targets"], bok["share_accounted"], bok["split_half"], 100 * bok["split_half_frac"], bok["k_eff"]))
                    lines.append("     mean arm by lam: " + "  ".join("%g:%.4f" % (l, A[:, i].mean()) for i, l in enumerate(lams)))
            g = {q: float(np.nanmean([r["%s@%g" % (q, out["lam_chosen_per_fold"][int(f)])] for r, f in zip(rr, folds)])) for q in ("gam_prob", "cos_prob", "gam_loc", "cos_loc", "mae")}
        else:
            for k in ("arm", "cloud", "sel"):
                vals[k] = np.array([r[k] for r in rr])
            g = {q: float(np.nanmean([r[q] for r in rr])) for q in ("gam_prob", "cos_prob", "gam_loc", "cos_loc", "mae")}
        out["stats"] = {}
        for k, lab in (("arm", "BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)"), ("cloud", "POINT CLOUD (3.0483 basis)"), ("sel", "SELECTION argmin K=500 (3.4540 basis)")):
            b = np.array([r[k] for r in base])
            c = ST.compare(vals[k], b, folds, names=pdbs, label="%s vs shipped -- %s" % (rung, lab))
            out["stats"][k] = c; lines.append(ST.fmt(c))
        out["progress"] = g
        lines.append("  gamma-equivalent: gam_eff prob-space %+.4f at cos %+.3f ; loc-space %+.4f at cos %+.3f ; MAE %.4f (diagnostic only)."
                     % (g["gam_prob"], g["cos_prob"], g["gam_loc"], g["cos_loc"], g["mae"]))
        lines.append("  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.")
        lines.append("  folds same sign (arm): %d/5 ; verdict (arm): %s" % (out["stats"]["arm"]["folds_same_sign"], out["stats"]["arm"]["verdict"]))
    text = "\n".join(lines)
    print(text)
    out["text"] = text
    ST.save_atomic(os.path.join(HERE, "results", "p_ladder_report_%s_s%d.json" % (rung, seed)), out, module_file=__file__)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rung", required=True); ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(); report(a.rung, a.seed)
