"""s17/sel_cal.py -- TARGET-LEVEL CALIBRATION, and the control that decides whether it exists.

Sprint 17 section 14, PREREG E4.  Two independent lines say the missing quantity is a
per-target scalar supplied at inference: in-band ordering is 0.986 within a target and 0.600
across, and Sprint 16's compactness residual has real magnitude with a per-target SIGN that
four native-free estimators predicted worse than a constant guess (0.357-0.413 against 0.516).

THE TRAP THIS MODULE EXISTS TO AVOID, and it caught this workstream inside an hour.

Given V selectors and T targets, the quantity `mean_t min_v RMSD(t, v)` looks like the ceiling
of per-target calibration.  It is not.  The minimum of V noisy draws falls with V whether or
not any per-target structure exists.  On the 58-variant sweep that statistic reads 2.482 A
against a shipped 3.454 -- an apparent 0.972 A of headroom.  **The matched zero-information
reference -- the minimum over V INDEPENDENT RANDOM selections from the same candidate set --
reads 2.148 A, which is BETTER.**  The 58 real variants are correlated enough with one another
that their per-target minimum does not even reach the noise floor.  There is no headroom
there at all; there is a min-of-N artefact with the sign of a discovery.

So every per-target-choice ceiling in this module is printed beside `minN_null`, the same
statistic over the same number of zero-information selectors, and a ceiling that does not
clear its null is reported as absent rather than as a target.

WHAT IS ACTUALLY TESTED.  A small, pre-declared arm set (below), and three baselines the
calibrator must beat, exactly as pre-registered:

    uncalibrated   the single arm the pipeline would use anyway
    constant       the single best arm chosen leave-fold-out -- a ZERO-INFORMATION calibration
    linear         a linear model on the same features, so a nonlinear one must earn itself

FEATURES ARE STRICTLY NATIVE-FREE and are listed in `FEATURES` so an auditor can check every
one of them without reading the code.  Labels (native RMSD) enter only inside training folds.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from s17 import sel_lib as L             # noqa: E402
from s17 import sel_obj as O             # noqa: E402

HYD = set("AVILMFWC")
CHG = set("DEKR")
POL = set("STNQYH")

#: every feature, in words, so the leakage claim is checkable without reading code
FEATURES = [
    "n                chain length",
    "logU             log universe size",
    "f_hydrophobic    fraction of AVILMFWC",
    "f_charged        fraction of DEKR",
    "f_polar          fraction of STNQYH",
    "f_gly            fraction of G",
    "f_pro            fraction of P",
    "net_charge       (K+R-D-E)/n",
    "hyd_moment       |sum exp(2*pi*i*k/3.6) * hydropathy_k| / n   (helical patterning)",
    "hyd_run          longest run of consecutive hydrophobics / n",
    "pred_helix       distogram-predicted helix content, from E[d(i,i+4)]",
    "pred_strand      distogram-predicted strand content, from E[d(i,i+4)]",
    "pred_rg          distogram-predicted radius of gyration, closed form",
    "pred_rg_z        the same, residualised on length",
    "sd_mean          mean predicted pair sd  (the model's own confidence)",
    "sd_max           max predicted pair sd",
    "ent_mean         mean predicted-distribution entropy",
    "sim_mean         mean BLOSUM similarity of the retrieved pool",
    "sim_gap          (top sim - median sim) / sd sim",
    "f_org            fraction of the pool from the peptide database, not protein fragments",
    "sc_mean          mean shipped score over the candidate set",
    "sc_sd            sd of the shipped score",
    "sc_skew          skew of the shipped score",
    "sc_gap           (median - min) / sd of the shipped score  -- best-vs-median gap",
    "pool_disp        mean typicality: mean_p mean_k |D_kp - setmean_p|",
    "pool_spread      mean pairwise CA-RMSD within a random 40-member subsample",
    "agree_top        mean pairwise CA-RMSD within the score's own top 25 (candidate agreement)",
]

ARMS = ("dist", "cons25", "cons75", "cons150", "blend15")


def _hydropathy():
    kd = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
              L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
    return kd


def target_features(p, sc, Wk, u):
    """All native-free per-target features.  `p` is a sel_lib pack, `sc` the shipped score."""
    seq, n = p["seq"], p["n"]
    kd = _hydropathy()
    f = {}
    f["n"] = float(n)
    f["logU"] = float(np.log(p["U"]))
    f["f_hydrophobic"] = sum(c in HYD for c in seq) / n
    f["f_charged"] = sum(c in CHG for c in seq) / n
    f["f_polar"] = sum(c in POL for c in seq) / n
    f["f_gly"] = seq.count("G") / n
    f["f_pro"] = seq.count("P") / n
    f["net_charge"] = (seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")) / n
    h = np.array([kd.get(c, 0.0) for c in seq])
    ang = 2 * np.pi * np.arange(n) / 3.6
    f["hyd_moment"] = float(abs((h * np.exp(1j * ang)).sum()) / n)
    run = best = 0
    for c in seq:
        run = run + 1 if c in HYD else 0
        best = max(best, run)
    f["hyd_run"] = best / n

    sep = p["sep"]; exp = p["expected"]
    m4 = sep == 4
    d4 = exp[m4] if m4.any() else np.array([10.0])
    f["pred_helix"] = float(np.mean(np.abs(d4 - 6.2) < 1.5))
    f["pred_strand"] = float(np.mean(d4 > 11.0))
    # Rg^2 = (1/(2 N^2)) sum_ij d_ij^2 ; |i-j|<2 terms filled with ideal CA geometry
    s2 = float((exp ** 2).sum() * 2) + 2 * (n - 1) * 3.8 ** 2
    f["pred_rg"] = float(np.sqrt(max(s2 / (2 * n * n), 1e-6)))
    f["pred_rg_z"] = f["pred_rg"] / (n ** 0.4)
    f["sd_mean"] = float(p["sd"].mean()); f["sd_max"] = float(p["sd"].max())
    ent = -(p["prob"] * np.log(np.maximum(p["prob"], 1e-12))).sum(1)
    f["ent_mean"] = float(ent.mean())

    sim = np.asarray(u["sim"], float)[np.asarray(u["order"], int)][:p["k"]]
    f["sim_mean"] = float(sim.mean())
    f["sim_gap"] = float((sim[0] - np.median(sim)) / max(sim.std(), 1e-6))
    f["f_org"] = float(np.asarray(u["org"], bool)[np.asarray(u["order"], int)][:p["k"]].mean())

    f["sc_mean"] = float(sc.mean()); f["sc_sd"] = float(sc.std())
    z = (sc - sc.mean()) / max(sc.std(), 1e-9)
    f["sc_skew"] = float((z ** 3).mean())
    f["sc_gap"] = float((np.median(sc) - sc.min()) / max(sc.std(), 1e-9))

    D = np.asarray(p["D"], float)
    f["pool_disp"] = float(np.abs(D - D.mean(0)[None, :]).mean())
    rng = SD.stable_rng(p["pdb"], "s17feat")
    sub = rng.choice(len(Wk), min(40, len(Wk)), replace=False)
    f["pool_spread"] = float(I.pairwise_rmsd(Wk[sub]).mean())
    top = np.argsort(sc, kind="stable")[:25]
    f["agree_top"] = float(I.pairwise_rmsd(Wk[top]).mean())
    return f


def _cons(W, sc, m):
    top = np.argsort(sc, kind="stable")[:min(m, len(sc))]
    P = I.pairwise_rmsd(W[top])
    return int(top[int(np.argmin(P.mean(1)))])


def run(K=500, targets=None, verbose=True):
    tg = targets if targets is not None else L.targets()
    rows = []
    for q, t in enumerate(tg):
        p = L.pack(t["pdb"], K, want=("D", "W"))
        u = I.load_univ(t["pdb"])
        sc = L.shipped(p)
        rr = p["rr"]; W = p["W"]
        c = O.context(p)
        bl = O.blend_scores(c, betas=(0.15,), scale="pred")["blendpred_0.15"]
        arms = {
            "dist": L.sel_of(sc, rr),
            "cons25": float(rr[_cons(W, sc, 25)]),
            "cons75": float(rr[_cons(W, sc, 75)]),
            "cons150": float(rr[_cons(W, sc, 150)]),
            "blend15": L.sel_of(bl, rr),
        }
        rows.append({"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]),
                     "k": int(p["k"]), "oracle": float(rr.min()), "random": float(rr.mean()),
                     "arms": arms, "feat": target_features(p, sc, W, u)})
        if verbose and (q + 1) % 25 == 0:
            print(f"  {q+1}/{len(tg)}", flush=True)
            json.dump({"rows": rows, "features": FEATURES},
                      open(os.path.join(L.RESULTS, "sel_cal.json"), "w"))
    json.dump({"rows": rows, "features": FEATURES},
              open(os.path.join(L.RESULTS, "sel_cal.json"), "w"))
    return rows


# --------------------------------------------------------------------------- models
def _ridge(X, y, lam):
    X = np.hstack([np.ones((len(X), 1)), X])
    A = X.T @ X + lam * np.eye(X.shape[1])
    A[0, 0] -= lam
    return np.linalg.solve(A, X.T @ y)


def _apply(w, X):
    return np.hstack([np.ones((len(X), 1)), X]) @ w


def minN_null(rows, V, tag):
    """min over V zero-information selections from the same candidate set."""
    out = []
    for r in rows:
        u = I.load_univ(r["pdb"])
        o = np.asarray(u["order"], int)
        rr = np.asarray(u["rr"], float)[o][:r["k"]]
        rng = SD.stable_rng(r["pdb"], "s17minN", tag)
        out.append(np.mean([rr[rng.integers(0, len(rr), V)].min() for _ in range(64)]))
    return np.array(out)


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(L.RESULTS, "sel_cal.json")))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    fk = sorted(rows[0]["feat"].keys())
    X = np.array([[r["feat"][k] for k in fk] for r in rows], float)
    X = (X - X.mean(0)) / np.maximum(X.std(0), 1e-9)
    A = np.array([[r["arms"][a] for a in ARMS] for r in rows], float)
    oracle = np.array([r["oracle"] for r in rows])
    rand = np.array([r["random"] for r in rows])
    dist = A[:, ARMS.index("dist")]
    n = len(rows)

    print(f"\n{'='*100}\nTARGET-LEVEL CALIBRATION   n = {n} targets, K = {rows[0]['k']}, "
          f"{len(ARMS)} pre-declared arms\n{'='*100}")
    print("ORACLE CEILING (best member)          %7.3f" % oracle.mean())
    print("MATCHED RANDOM (one draw)             %7.3f" % rand.mean())
    for a in ARMS:
        v = A[:, ARMS.index(a)]
        r = L.report_pair(f"{a} vs dist", v, dist, fold)
        print(f"  {a:<10}{v.mean():>8.3f}   vs dist {r['diff']:+.3f} "
              f"[{r['lo']:+.3f},{r['hi']:+.3f}]  {r['W']}W/{r['L']}L")

    print("\nA. THE PER-TARGET-CHOICE CEILING, WITH ITS ZERO-INFORMATION NULL")
    nul = minN_null(rows, len(ARMS), "arms")
    print(f"  per-target ORACLE over the {len(ARMS)} arms   {A.min(1).mean():7.3f}   [ORACLE]")
    print(f"  min over {len(ARMS)} ZERO-INFO selections      {nul.mean():7.3f}   [zero-information]")
    print(f"  the ceiling clears its null by       {nul.mean()-A.min(1).mean():+7.3f} A "
          f"({'REAL headroom' if A.min(1).mean() < nul.mean() else 'NO headroom -- a min-of-N artefact'})")

    print("\nB. THE THREE BASELINES THE CALIBRATOR MUST BEAT  (PREREG E4)")
    unc = dist
    lfo_const, choice, _ = L.lfo_variant(A, fold)
    print(f"  uncalibrated (dist)                  {unc.mean():7.3f}")
    print(f"  constant   (best arm, leave-fold-out) {lfo_const.mean():7.3f}   "
          f"choices {dict((k, ARMS[v]) for k, v in choice.items())}")

    # linear / calibrated: leave-fold-out ridge predicting each arm's margin over dist
    for lam in (1.0, 10.0, 100.0):
        pred_sel = np.empty(n)
        for g in np.unique(fold):
            tr = fold != g
            M = np.empty((int((~tr).sum()), len(ARMS)))
            for a in range(len(ARMS)):
                w = _ridge(X[tr], A[tr, a] - dist[tr], lam)
                M[:, a] = _apply(w, X[~tr])
            pred_sel[~tr] = A[~tr, np.argmin(M, axis=1)]
        r1 = L.report_pair("calibrated vs uncalibrated", pred_sel, unc, fold)
        r2 = L.report_pair("calibrated vs constant", pred_sel, lfo_const, fold)
        print(f"  linear ridge lam={lam:<6g}            {pred_sel.mean():7.3f}")
        print("     " + L.fmt_pair(r1))
        print("     " + L.fmt_pair(r2))

    print("\nC. DO ANY FEATURES CARRY THE PER-TARGET MARGIN AT ALL?")
    print("   Spearman of each feature against (cons75 - dist), the per-target margin whose")
    print("   SIGN a calibrator would have to supply.  Null band at n=126 is about +-0.175.")
    marg = A[:, ARMS.index("cons75")] - dist
    rs = [(O.spearman(X[:, q], marg), fk[q]) for q in range(len(fk))]
    for v, k in sorted(rs, key=lambda z: -abs(z[0]))[:10]:
        print(f"    {k:<16}{v:+.3f}")
    print(f"    sign of the margin: {int((marg<0).sum())} favour cons75, "
          f"{int((marg>0).sum())} favour dist, {int((marg==0).sum())} tied")
    print(f"    constant-sign baseline (always the majority arm): "
          f"{max((marg<0).mean(), (marg>0).mean()):.3f}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
        report()
