"""s26/a_strain_vs_spread.py -- lane A: adversary check of L53 (strain difficulty).

L53 reports partial Spearman(+0.433) between `amber_moved` (native-free) and `rmsd_arm` (ORACLE
label), partialled on n and Rg. The prereg carried no pool-disagreement control. The obvious
native-free difficulty proxy is the shipped top-75's own pairwise CA-RMSD spread (the members'
disagreement, no native). This script asks: does `moved` carry signal BEYOND that spread, or is
it the spread in disguise? ORACLE DIAGNOSTIC: `rmsd_arm` is read as the label only; nothing is
selected or tuned. Reads the production cache and the universes; gate is open (L33).

    python s26/jobrun.py --agent A --tag CPU --name a_strain_vs_spread --est-ram 0.4 -- \
        python s26/a_strain_vs_spread.py
"""
from __future__ import annotations
import json, os, sys, glob
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from scipy.stats import spearmanr                    # noqa: E402
OUT = os.path.join(ROOT, "s26", "results", "a_strain_vs_spread.json")


def residualise(y, X):
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ beta


def rho_ci(x, y, folds, rng, nboot=2000):
    n = len(x)
    r = float(spearmanr(x, y)[0])
    bi = np.array([spearmanr(x[i], y[i])[0] for i in (rng.integers(0, n, n) for _ in range(nboot))])
    F = np.array(sorted(set(folds.tolist())))
    bf = []
    for _ in range(nboot):
        pick = np.concatenate([np.flatnonzero(folds == q) for q in rng.choice(F, len(F), replace=True)])
        bf.append(spearmanr(x[pick], y[pick])[0])
    bf = np.array(bf)
    perm = np.array([spearmanr(x, y[rng.permutation(n)])[0] for _ in range(nboot)])
    per_fold = {int(q): float(spearmanr(x[folds == q], y[folds == q])[0]) for q in F}
    return {"rho": r, "ci_iid": [float(np.percentile(bi, 2.5)), float(np.percentile(bi, 97.5))],
            "ci_fold": [float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))],
            "perm_p": float((np.abs(perm) >= abs(r)).mean()), "per_fold": per_fold,
            "folds_same_sign": int(sum(np.sign(v) == np.sign(r) for v in per_fold.values()))}


def main():
    rng = np.random.default_rng(2653)
    rows = []
    for t in I.targets():
        pdb = t["pdb"]
        rec = I.shipped_record(pdb)
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        sub = np.asarray(rec["sub"], int)
        top = W[sub]
        Pm = I.pairwise_rmsd(top)                       # native-free: members against members
        iu = np.triu_indices(len(sub), 1)
        ca = np.asarray(rec["ca"], float)
        rg = float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean()))
        rows.append({"pdb": pdb, "fold": int(t["fold"]), "n": int(t["n"]), "rg": rg,
                     "spread_mean": float(Pm[iu].mean()), "spread_median": float(np.median(Pm[iu])),
                     "medoid_min": float((Pm.sum(1) / (len(sub) - 1)).min()),
                     "moved": float(rec["amber_moved"]), "log_e0": float(np.log10(max(rec["amber_e0"], 1.0))),
                     "rmsd_arm": float(rec["rmsd_arm"])})
        print(f"{pdb} spread {rows[-1]['spread_mean']:.3f} moved {rows[-1]['moved']:.3f} arm {rows[-1]['rmsd_arm']:.3f}", flush=True)
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    y = np.array([r["rmsd_arm"] for r in rows])
    X = np.column_stack([[r["n"] for r in rows], [r["rg"] for r in rows]])
    out = {"n": len(rows), "label": "ORACLE DIAGNOSTIC: native-free signals against rmsd_arm as label",
           "raw": {}, "partial_n_rg": {}, "partial_n_rg_spread": {}}
    sig = {k: np.array([r[k] for r in rows]) for k in ("spread_mean", "spread_median", "medoid_min", "moved", "log_e0")}
    y_res = residualise(y, X)
    for k, v in sig.items():
        out["raw"][k] = rho_ci(v, y, folds, rng)
        out["partial_n_rg"][k] = rho_ci(residualise(v, X), y_res, folds, rng)
    # moved and log_e0 beyond the spread
    Xs = np.column_stack([X, sig["spread_mean"]])
    y_res2 = residualise(y, Xs)
    for k in ("moved", "log_e0"):
        out["partial_n_rg_spread"][k] = rho_ci(residualise(sig[k], Xs), y_res2, folds, rng)
    out["rho_moved_spread"] = float(spearmanr(sig["moved"], sig["spread_mean"])[0])
    for blk in ("raw", "partial_n_rg", "partial_n_rg_spread"):
        for k, q in out[blk].items():
            print(f"{blk:22s} {k:14s} rho {q['rho']:+.3f}  iid {np.round(q['ci_iid'], 3).tolist()}  fold {np.round(q['ci_fold'], 3).tolist()}  perm p {q['perm_p']:.4f}  folds {q['folds_same_sign']}/5")
    print("rho(moved, spread_mean) =", round(out["rho_moved_spread"], 3))
    ST.save_atomic(OUT, {"summary": out, "rows": rows}, complete_keys=("pdb", "spread_mean", "moved", "rmsd_arm"),
                   rows=rows, n_expected=126, module_file=__file__)
    print("wrote", os.path.relpath(OUT, ROOT))


if __name__ == "__main__":
    main()
