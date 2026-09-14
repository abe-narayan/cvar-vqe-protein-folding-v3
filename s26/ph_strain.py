"""s26/ph_strain.py -- THE RELAXATION'S OWN STRAIN AS A NATIVE-FREE DIFFICULTY SIGNAL (lane PH).

Pre-registered in `s26/PREREG_strain_difficulty.md`. Gated: reads `rmsd_arm` as the ORACLE label
of the already-emitted built chain. Reads only the 126 production records. Under a minute.

    python s26/ph_strain.py [--rep]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s26 import ph_lib as L                                   # noqa: E402

OUT = os.path.join(L.RESULTS, "ph_strain.json")
NBOOT = 4000
SIGNALS = ("log_e0", "log_drop", "moved", "strain_after")


def spearman(a, b) -> float:
    from scipy.stats import spearmanr
    return float(spearmanr(a, b)[0])


def residualise(y, X):
    """y minus its least-squares fit on the columns of X (with intercept)."""
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ beta


def rho_with_cis(x, y, folds, rng, label):
    """Spearman with iid and fold-clustered bootstrap CIs, a permutation p and per-fold signs."""
    n = len(x)
    r = spearman(x, y)
    bi = np.array([spearman(x[i], y[i]) for i in (rng.integers(0, n, n) for _ in range(NBOOT))])
    F = np.array(sorted(set(folds.tolist())))
    bf = []
    for _ in range(NBOOT):
        pick = np.concatenate([np.flatnonzero(folds == q) for q in rng.choice(F, len(F), replace=True)])
        bf.append(spearman(x[pick], y[pick]))
    bf = np.array(bf)
    perm = np.array([spearman(x, y[rng.permutation(n)]) for _ in range(NBOOT)])
    p = float((np.abs(perm) >= abs(r)).mean())
    per_fold = {int(q): spearman(x[folds == q], y[folds == q]) for q in F}
    same = int(sum(np.sign(v) == np.sign(r) for v in per_fold.values()))
    return {"label": label, "rho": r, "ci95_iid": [float(np.percentile(bi, 2.5)), float(np.percentile(bi, 97.5))],
            "ci95_fold": [float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))],
            "perm_p": p, "per_fold": per_fold, "folds_same_sign": same, "n": int(n)}


def fisher_top_quartile(signal, is_fail, pdbs):
    from scipy.stats import fisher_exact
    k = len(signal) // 4
    top = set(np.argsort(-signal, kind="stable")[:k].tolist())
    a = sum(1 for i in top if is_fail[i]); b = k - a
    c = int(is_fail.sum()) - a; d = len(signal) - k - c
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {"top_k": k, "fail_in_top": a, "fail_total": int(is_fail.sum()), "odds": float(odds),
            "p_one_sided": float(p), "top_targets": sorted(pdbs[i] for i in top)}


def main(rep: bool = False) -> dict:
    L.require_gate("ph_strain")
    from s12 import instrument as I
    tg = L.targets()
    if rep:
        tg = list(reversed(tg))
    rows = []
    for t in tg:
        rec = L.prod_record_oracle(t["pdb"])
        ca = np.asarray(rec["ca"], float)
        rows.append({"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]),
                     "rmsd_arm": float(rec["rmsd_arm"]), "e0": float(rec["amber_e0"]), "e1": float(rec["amber_e1"]),
                     "moved": float(rec["amber_moved"]), "strain_after": float(rec["amber_strain_after"]),
                     "rg": float(L.rg_of(ca)), "fail18": t["pdb"] in I.FAIL18})
    rows = sorted(rows, key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]
    folds = L.folds_of(pdbs)
    y = np.array([r["rmsd_arm"] for r in rows])
    sig = {"log_e0": np.log10(np.maximum([r["e0"] for r in rows], 1.0)),
           "log_drop": np.log10(np.maximum([r["e0"] - r["e1"] for r in rows], 1.0)),
           "moved": np.array([r["moved"] for r in rows]),
           "strain_after": np.array([r["strain_after"] for r in rows])}
    X = np.column_stack([[r["n"] for r in rows], [r["rg"] for r in rows]]).astype(float)
    y_res = residualise(y, X)
    rng = L.stable_rng("strain", "rep" if rep else "main")
    out = {"n": len(rows), "rep": rep, "raw": {}, "partial_n_rg": {}, "fisher_fail18": {},
           "confounds": {"rho_n_rmsd": spearman(X[:, 0], y), "rho_rg_rmsd": spearman(X[:, 1], y),
                         "rho_n_log_e0": spearman(X[:, 0], sig["log_e0"]),
                         "rho_rg_log_e0": spearman(X[:, 1], sig["log_e0"])}}
    is_fail = np.array([r["fail18"] for r in rows])
    for k in SIGNALS:
        out["raw"][k] = rho_with_cis(sig[k], y, folds, rng, f"raw rho({k}, rmsd_arm)")
        out["partial_n_rg"][k] = rho_with_cis(residualise(sig[k], X), y_res, folds, rng,
                                              f"partial rho({k}, rmsd_arm | n, Rg)")
        out["fisher_fail18"][k] = fisher_top_quartile(sig[k], is_fail, pdbs)
    # the falsifier, on the partial correlations, Bonferroni over four signals
    hits = []
    for k in SIGNALS:
        q = out["partial_n_rg"][k]
        ok = (abs(q["rho"]) >= 0.25 and (q["ci95_fold"][0] > 0 or q["ci95_fold"][1] < 0)
              and q["folds_same_sign"] == 5 and q["perm_p"] < 0.05 / 4)
        if ok:
            hits.append(k)
    out["falsifier_clears"] = bool(hits); out["signals_clearing"] = hits
    for kind in ("raw", "partial_n_rg"):
        for k in SIGNALS:
            q = out[kind][k]
            print(f"{q['label']:<45} rho {q['rho']:+.3f}  iid [{q['ci95_iid'][0]:+.3f},{q['ci95_iid'][1]:+.3f}]"
                  f"  fold [{q['ci95_fold'][0]:+.3f},{q['ci95_fold'][1]:+.3f}]  perm p {q['perm_p']:.4f}"
                  f"  folds same sign {q['folds_same_sign']}/5")
    for k in SIGNALS:
        f = out["fisher_fail18"][k]
        print(f"FAIL18 in top quartile of {k:<13}: {f['fail_in_top']}/{f['top_k']} (of {f['fail_total']})  odds {f['odds']:.2f}  p {f['p_one_sided']:.3f}")
    print("confounds", {k: round(v, 3) for k, v in out["confounds"].items()})
    print("FALSIFIER CLEARS:", out["falsifier_clears"], hits)
    L.save(OUT.replace(".json", "_rep.json") if rep else OUT,
           {"what": "relaxation strain scalars vs the built chain's ORACLE error; calibration only", "rows": rows,
            "summary": out}, rows=rows, complete_keys=("pdb", "rmsd_arm", "e0", "e1", "moved", "strain_after"),
           n_expected=126, module_file=__file__)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rep", action="store_true")
    a = ap.parse_args()
    main(rep=a.rep)
