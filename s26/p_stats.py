"""s26/p_stats.py -- lane P's shared nested leave-fold-out ridge machinery (B3, C4, C5).

Dependency-free (numpy only) so the same closed-form ridge is used in every lane-P file:

    ridge_fit(X, y, alpha)            -> (w, b, mu, sd): standardise on the TRAINING rows only
    nested_predict(X, Y, folds, ...)  -> held-out predictions for every target, alpha chosen by
                                         inner leave-one-fold-out over the training folds
    balanced_accuracy, permutation_null

Every router in the record that used a tree ensemble overfit (S22 L7: in-fold -0.276 -> held-out
+0.079); S22 L10's bound says nothing richer than a linear router is learnable at n = 126. So
this is deliberately linear, with the alpha grid and the inner CV declared once, here.
"""
from __future__ import annotations

import numpy as np

ALPHAS = np.logspace(-2, 4, 25)


def ridge_fit(X, Y, alpha):
    """Closed-form ridge with intercept on standardised columns.  Y may be (n,) or (n, k)."""
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    mu = X.mean(0); sd = X.std(0) + 1e-9
    Z = (X - mu) / sd
    ym = Y.mean(0)
    A = Z.T @ Z + alpha * np.eye(Z.shape[1])
    W = np.linalg.solve(A, Z.T @ (Y - ym))
    return {"W": W, "b": ym, "mu": mu, "sd": sd}


def ridge_predict(m, X):
    return ((np.asarray(X, float) - m["mu"]) / m["sd"]) @ m["W"] + m["b"]


def _score(pred, Y, kind):
    if kind == "class":
        return balanced_accuracy(np.sign(Y), np.sign(pred))
    return -float(((pred - Y) ** 2).mean())


def nested_predict(X, Y, folds, kind="reg", alphas=ALPHAS):
    """Held-out prediction for every row.  For each outer fold f: alpha is chosen by inner
    leave-one-fold-out over the OTHER folds (maximising balanced accuracy for `class`, minimising
    MSE for `reg`), the model is refit on all other folds at that alpha, and fold f is predicted.
    Returns (pred, chosen_alpha_per_fold)."""
    X = np.asarray(X, float); Y = np.asarray(Y, float); folds = np.asarray(folds)
    F = sorted(set(folds.tolist()))
    pred = np.zeros(Y.shape, float); chosen = {}
    for f in F:
        tr = folds != f
        inner = [g for g in F if g != f]
        best, best_a = -np.inf, alphas[0]
        for a in alphas:
            s = []
            for g in inner:
                itr = tr & (folds != g); ite = folds == g
                m = ridge_fit(X[itr], Y[itr], a)
                s.append(_score(ridge_predict(m, X[ite]), Y[ite], kind))
            s = float(np.mean(s))
            if s > best + 1e-12:
                best, best_a = s, a
        m = ridge_fit(X[tr], Y[tr], best_a)
        pred[folds == f] = ridge_predict(m, X[folds == f])
        chosen[int(f)] = float(best_a)
    return pred, chosen


def balanced_accuracy(y, yhat):
    y = np.sign(np.asarray(y, float)); yhat = np.sign(np.asarray(yhat, float))
    yhat[yhat == 0] = 1.0
    out = []
    for c in (-1.0, 1.0):
        m = y == c
        if m.any():
            out.append(float((yhat[m] == c).mean()))
    return float(np.mean(out)) if out else float("nan")


def permutation_null(X, Y, folds, kind, n_perm, rng, alphas=ALPHAS):
    """The same nested procedure on labels permuted across targets."""
    vals = []
    for _ in range(n_perm):
        Yp = Y[rng.permutation(len(Y))]
        p, _ = nested_predict(X, Yp, folds, kind, alphas)
        vals.append(balanced_accuracy(Yp, p) if kind == "class" else -float(((p - Yp) ** 2).mean()))
    return np.array(vals)


def selftest():
    rng = np.random.default_rng(0)
    n, p = 126, 8
    folds = np.repeat(np.arange(5), 26)[:n]
    X = rng.normal(size=(n, p))
    w = rng.normal(size=p)
    y_sig = np.sign(X @ w + 0.5 * rng.normal(size=n))
    y_rnd = np.sign(rng.normal(size=n))
    ps, _ = nested_predict(X, y_sig, folds, "class")
    pr, _ = nested_predict(X, y_rnd, folds, "class")
    acc_s, acc_r = balanced_accuracy(y_sig, ps), balanced_accuracy(y_rnd, pr)
    null = permutation_null(X, y_rnd, folds, "class", 20, rng)
    print("  planted signal: held-out balanced acc %.3f ; random labels %.3f ; null 95th pct %.3f"
          % (acc_s, acc_r, np.percentile(null, 95)))
    assert acc_s > np.percentile(null, 95) and acc_s > 0.75
    assert acc_r <= max(np.percentile(null, 99), 0.62)
    d = X @ w + 0.3 * rng.normal(size=n)
    pd, _ = nested_predict(X, d, folds, "reg")
    r2 = 1 - ((d - pd) ** 2).sum() / ((d - d.mean()) ** 2).sum()
    print("  planted regression: held-out R2 %.3f" % r2)
    assert r2 > 0.8
    print("  p_stats selftest OK")


if __name__ == "__main__":
    selftest()
