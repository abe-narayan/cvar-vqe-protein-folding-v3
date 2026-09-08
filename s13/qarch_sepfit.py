"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 2c -- Legacy vs AMBER separability, matched.

The exact ANOVA of section 2 needs the full 4^9 factorial, which AMBER cannot afford.  This
file gets the same quantity -- the ORDER-1 (separable) variance share -- from the AMBER
single points ALREADY CACHED by `qarch_enum` on its 1,194-configuration uniform subsample,
with no new energy evaluations at all.

Estimator.  The best additive approximation of F under the uniform product measure is the
least-squares fit of one-hot residue-state indicators, so its (cross-validated) R^2 estimates
the ANOVA's V1 exactly.  Order 2 adds the pair indicators.

VALIDATION IS THE POINT: Legacy has BOTH the exact enumeration V1 and the sampled estimate on
the identical 1,194 configurations, so the estimator is calibrated against ground truth before
it is trusted on AMBER.  AMBER is rank-transformed first, because its raw single points span
20 orders of magnitude (section 3d) and a least-squares fit on them is meaningless; Legacy is
reported both raw and rank-transformed so the comparison is like for like.

Output: `s13/results/qarch_sepfit.json`.

    python -m s13.qarch_sepfit
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402


def design(S, k, order=1):
    """One-hot (and pair) indicator design matrix for state configurations S (B, n)."""
    S = np.asarray(S, int)
    B, n = S.shape
    cols = [np.ones((B, 1))]
    for i in range(n):
        M = np.zeros((B, k - 1))
        for a in range(k - 1):
            M[:, a] = (S[:, i] == a)
        cols.append(M)
    if order >= 2:
        for i, j in itertools.combinations(range(n), 2):
            M = np.zeros((B, (k - 1) ** 2))
            c = 0
            for a in range(k - 1):
                for b in range(k - 1):
                    M[:, c] = (S[:, i] == a) & (S[:, j] == b)
                    c += 1
            cols.append(M)
    return np.hstack(cols)


def cv_r2(X, y, folds=5, seed=0, ridge=1e-8):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(y))
    parts = np.array_split(idx, folds)
    pred = np.empty(len(y))
    for f in range(folds):
        te = parts[f]
        tr = np.concatenate([parts[g] for g in range(folds) if g != f])
        A = X[tr]
        beta = np.linalg.solve(A.T @ A + ridge * np.eye(X.shape[1]), A.T @ y[tr])
        pred[te] = X[te] @ beta
    ss = ((y - y.mean()) ** 2).sum()
    return float(1.0 - ((y - pred) ** 2).sum() / ss)


def main():
    pdbs = sorted(f[len("qarch_enum_"):-4] for f in os.listdir(Q.RESULTS)
                  if f.startswith("qarch_enum_") and f.endswith(".npz"))
    rows = []
    for p in pdbs:
        z = np.load(os.path.join(Q.RESULTS, f"qarch_enum_{p}.npz"))
        n, k = int(z["n"]), int(z["k"])
        a_idx = np.asarray(z["amber_idx"], int)
        m = np.asarray(z["amber_kind"], int) == 0
        idx = a_idx[m]
        S = np.array(np.unravel_index(idx, (k,) * n)).T
        X1 = design(S, k, 1); X2 = design(S, k, 2)
        A = np.asarray(z["amber_total"], float)[m]
        L = np.asarray(z["legacy"], float)[idx]
        R = np.asarray(z["rmsd"], float)[idx]
        rec = {"pdb": p, "n": n, "k": k, "n_samples": int(len(idx)),
               "n_params_order1": int(X1.shape[1]), "n_params_order2": int(X2.shape[1])}
        for nm, y in (("legacy_raw", L), ("legacy_rank", Q._rank(L)),
                      ("amber_rank", Q._rank(A)), ("ORACLE_rmsd", R)):
            rec[nm] = {"cv_r2_order1": cv_r2(X1, np.asarray(y, float)),
                       "cv_r2_order2": cv_r2(X2, np.asarray(y, float))}
        rows.append(rec)
        print(f"  {p}: legacy_rank V1~{rec['legacy_rank']['cv_r2_order1']:+.3f} "
              f"amber_rank V1~{rec['amber_rank']['cv_r2_order1']:+.3f} "
              f"(n={len(idx)}, p1={X1.shape[1]}, p2={X2.shape[1]})", flush=True)
    Q.write("qarch_sepfit", {"what": "sampled separability of Legacy vs AMBER on identical "
                                     "cached configurations", "rows": rows})
    print(f"\n{'quantity':<16s} {'cv R2 order1':>13s} {'cv R2 order2':>13s}")
    for nm in ("legacy_raw", "legacy_rank", "amber_rank", "ORACLE_rmsd"):
        a = np.mean([r[nm]["cv_r2_order1"] for r in rows])
        b = np.mean([r[nm]["cv_r2_order2"] for r in rows])
        print(f"{nm:<16s} {a:13.3f} {b:13.3f}")
    return rows


if __name__ == "__main__":
    main()
