"""s22/errweight.py -- ERROR-PLACEMENT WEIGHTING, H_D* = sum_ij s_ij w_ij (dhat_ij - d_ij)^2
(PREREG_C.md Experiment 2, Priority 2 of s22/BRIEF.md).

Four PRE-REGISTERED, closed-form, native-free weighting schemes over pool candidate pairs
(i,j), i<j-1 (the distogram's own pair index, |i-j|>=2), all multiplying the existing
precision term w_ij = 1/sd_ij^2 (a Gaussian approximation to the shipped Bayes-risk score --
declared and soundness-gated below, since the shipped risk table has no closed pairwise-sum
decomposition to re-weight):

    s_control    1                    plain precision-weighted squared error (the null)
    s_local      1 / (j - i)          upweight LOCAL pairs (torsion-space locality theorem:
                                       d_ij depends on exactly j-i-1 torsions -- fewer torsions,
                                       more directly localisable per-torsion signal)
    s_global     (j - i)              upweight LONG-RANGE pairs (RMSD is a GLOBAL statistic;
                                       local pairs can be individually satisfied while the fold
                                       is globally wrong)
    s_discrim    Var_pool(d_ij)       upweight pairs where the retrieval POOL's own 500
                                       candidates DISAGREE -- native-free proxy for "this pair
                                       carries information that can tell candidates apart"

None is fit to RMSD; all are fixed closed-form functions of (i,j) or of the POOL's own realised
geometry (never the native).  No CV question therefore arises for this file.

Reports, per scheme, on n=126: the CERTIFIED OPTIMUM (pool argmin under H_D*) for mechanism, and
the DEPLOYED READOUT (top-75-by-H_D* coordinate average) for the primary -- never conflated, per
BRIEF SS7.  `s_control` is the null every structural scheme must beat, and is itself checked for
soundness against the shipped Bayes-risk score's own ranking (a Gaussian approximation is a new
functional fork beyond the project's existing squared-vs-Bayes one, so it is validated before any
s_ij != 1 result is trusted).

FALSIFIER (PREREG_C.md Sec 2): dead unless a structural scheme beats BOTH s_control AND the shipped
incumbent (3.048) on the DEPLOYED readout, past that comparison's own MDE, with a CI excluding zero.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

TOPM = 75
SCHEMES = ("s_control", "s_local", "s_global", "s_discrim")


def _sij(scheme, i, j, W=None):
    if scheme == "s_control":
        return np.ones(len(i))
    if scheme == "s_local":
        return 1.0 / (j - i).astype(float)
    if scheme == "s_global":
        return (j - i).astype(float)
    if scheme == "s_discrim":
        D = I.pair_dists(W, i, j)          # (K, npairs)
        return D.var(axis=0)
    raise ValueError(scheme)


def _save(obj):
    path = os.path.join(RESULTS, "errweight.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def run():
    tg = I.targets()
    rows = []
    print("targets: %d" % len(tg), flush=True)
    for c, t in enumerate(tg):
        pdb, n = t["pdb"], int(t["n"])
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(n)
        dhat = np.asarray(dg["expected"], float)
        sd = np.asarray(dg["sd"], float)
        sd2 = np.maximum(sd, 1e-3) ** 2
        D = I.pair_dists(W, i, j)              # (K, npairs)
        resid2 = (dhat[None, :] - D) ** 2
        rr = I.kabsch_rmsd_batch(W, nat)

        # shipped score, for the soundness check
        sc_shipped = np.asarray(I.shipped_score(dg, D), float)
        rank_shipped = np.argsort(sc_shipped, kind="stable")

        row = {"pdb": pdb, "n": n, "fold": int(t["fold"])}
        for scheme in SCHEMES:
            s_ij = _sij(scheme, i, j, W=W)
            weight = s_ij / sd2
            score = (resid2 * weight[None, :]).sum(1)
            order = np.argsort(score, kind="stable")
            row["%s_argmin" % scheme] = float(rr[order[0]])
            C, _b = I.coordinate_average(W[order[:TOPM]])
            row["%s_avg75" % scheme] = float(I.ca_rmsd(np.asarray(C, float), nat))
            if scheme == "s_control":
                # soundness check: rank correlation vs the shipped Bayes-risk score
                from scipy.stats import spearmanr
                row["control_vs_shipped_spearman"] = float(spearmanr(score, sc_shipped).statistic)
                row["control_argmin_eq_shipped_argmin"] = bool(order[0] == rank_shipped[0])
        rows.append(row)
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    ok = len(rows) == len(tg) and all(np.isfinite(r["s_control_avg75"]) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "schemes": list(SCHEMES)})
    report(rows)
    return rows


def _stat(d, rng, B=5000):
    d = np.asarray(d, float); k = len(d)
    se = float(d.std(ddof=1) / np.sqrt(k))
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return dict(mean=float(d.mean()), se=se, mde=2.8016 * se,
                lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)),
                w=int((d < 0).sum()), l=int((d > 0).sum()))


def _fmt(s):
    return "%+.4f SE %.4f MDE %.4f [%+.4f,%+.4f] %dW/%dL" % (s["mean"], s["se"], s["mde"], s["lo"], s["hi"], s["w"], s["l"])


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "errweight.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)     # noqa: E731
    rng = SD.stable_rng("s22", "errweight", "report")
    print("\nn = %d.  SOUNDNESS CHECK: s_control (Gaussian approx) vs shipped Bayes-risk score."
          % len(rows))
    print("  median Spearman(control_score, shipped_score) = %.4f   argmin agreement %.1f%%"
          % (np.median(g("control_vs_shipped_spearman")),
             100.0 * np.mean([r["control_argmin_eq_shipped_argmin"] for r in rows])))

    inc75 = None
    pg_path = os.path.join(RESULTS, "..", "..", "s21", "results", "poolgap.json")
    if os.path.exists(pg_path):
        pg = {r["pdb"]: r for r in json.load(open(pg_path))["rows"]}
        inc75 = np.array([pg[r["pdb"]]["avg_75"] for r in rows], float)
        print("  shipped incumbent (avg_75, cross-target-matched from poolgap.json): %.4f" % inc75.mean())

    print("\n  CERTIFIED OPTIMUM (pool argmin), each scheme vs s_control:")
    for scheme in SCHEMES:
        x = g("%s_argmin" % scheme)
        print("    %-12s mean %.4f  median %.4f" % (scheme, x.mean(), np.median(x)))
    ctrl_arg = g("s_control_argmin")
    for scheme in ("s_local", "s_global", "s_discrim"):
        x = g("%s_argmin" % scheme)
        s = _stat(x - ctrl_arg, rng)
        print("    %s - s_control:  %s" % (scheme, _fmt(s)))

    print("\n  DEPLOYED READOUT (top-75 average), each scheme vs s_control AND vs the shipped incumbent:")
    ctrl_avg = g("s_control_avg75")
    for scheme in SCHEMES:
        x = g("%s_avg75" % scheme)
        print("    %-12s mean %.4f  median %.4f" % (scheme, x.mean(), np.median(x)))
    for scheme in ("s_local", "s_global", "s_discrim"):
        x = g("%s_avg75" % scheme)
        s_ctrl = _stat(x - ctrl_avg, rng)
        print("    %-10s - s_control:            %s" % (scheme, _fmt(s_ctrl)))
        if inc75 is not None:
            s_inc = _stat(x - inc75, rng)
            print("    %-10s - shipped incumbent:    %s" % (scheme, _fmt(s_inc)))
    print("\n  FALSIFIER: dead unless a structural scheme beats BOTH s_control AND the shipped")
    print("  incumbent on the DEPLOYED readout past its own MDE with a CI excluding zero.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
