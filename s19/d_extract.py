"""AGENT D, Sprint 19 -- Block P6: is the mode information EXTRACTABLE by any native-free rule?

The coordinator asked for the one number that bounds the direction.  D's P2b measured
**-0.414 A [-0.502, -0.324]** of genuine mode information per multimodal pair (ORACLE nearest
mode minus its matched min-of-N null), while both native-free mode rules tried were WORSE than
the deployed mean.  This module asks whether ANY native-free function of the predicted
distribution's shape reaches it, and bounds the class with an ORACLE in-fold fit.

This is a PER-PAIR ERROR diagnostic, not an RMSD result.  BRIEF section 5 is explicit that a
lower aggregate error metric that does not lower Ca-RMSD is a negative result; P5 is the
structural test.  Nothing here is a predictive claim about RMSD.

ARMS -- every one is an aim point for the same pair, scored by |aim - d_true|

    mean_raw        the distogram's first moment, undebiased
    deployed        `mean - b(sep)`, the shipped leave-fold-out separation debias   <- THE BASELINE
    lfo_noshape     leave-fold-out ridge on sep + mean + sd only (no distribution shape)
    lfo_shape       leave-fold-out ridge on the full 17-bin vector + shape summaries
    lfo_gbm         leave-fold-out gradient boosting on the same features
    lfo_perm        the SAME model with labels permuted within the training fold -- zero-info
    ORACLE_infold   the same model fitted ON the held-out fold -- the CEILING of the rule class
    ORACLE_bestmode the nearest mode to the truth               } reproduced here so every
    minofN_null     matched-|offset| random-sign, minimised     } number sits on one scale

Fold assignment is the project's pinned 5 folds.  Seeding via `s15/seed.stable_rng`.

Also: block P2a robustness.  The cell-matched estimate was NOT MEASURED and the regression
adjustment was significant; the coordinator asked whether "at least as close" overstates it.
A non-parametric nearest-neighbour match on (separation exact, sd nearest, within target)
settles which reading the data supports.

Run:  python -m s19.d_extract
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I            # noqa: E402
from s15 import distcal as C               # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s19 import d_modality as M            # noqa: E402

OUT = os.path.join(HERE, "results", "D_P6_extract")
os.makedirs(OUT, exist_ok=True)
CENTRES, WIDTHS = M.CENTRES, M.WIDTHS


def build():
    tg = I.targets()
    data = C.gather(tg)                       # full target list -- the subset trap, section 7
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([data[p]["fold"] for p in pdbs], int)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = fn

    F, Y, T, FO, EXTRA = [], [], [], [], []
    for ti, p in enumerate(pdbs):
        d = data[p]
        z = np.load(os.path.join(ROOT, "s12", "cache", f"disto_{p}.npz"))
        prob = np.asarray(z["prob"], float)
        risk, grid = np.asarray(z["risk"], float), np.asarray(z["grid"], float)
        exp, sd, sep = d["dhat"], d["sd"], d["sep"]
        med = grid[np.argmin(risk, axis=1)]
        Mask = M.dens_modes(prob)
        Mp = M.prom_modes(prob)
        dens = np.where(Mask, prob / WIDTHS[None], -np.inf)
        o1 = np.argmax(dens, 1)
        d2 = dens.copy()
        d2[np.arange(len(o1)), o1] = -np.inf
        o2 = np.argmax(d2, 1)
        has2 = np.isfinite(d2[np.arange(len(o1)), o2])
        ent = -(prob * np.log(np.maximum(prob, 1e-12))).sum(1)
        mu3 = (prob * (CENTRES[None] - exp[:, None]) ** 3).sum(1) / np.maximum(sd, 1e-6) ** 3
        mu4 = (prob * (CENTRES[None] - exp[:, None]) ** 4).sum(1) / np.maximum(sd, 1e-6) ** 4
        m1a = (prob * (np.abs(CENTRES[None] - exp[:, None]) <= 1.0)).sum(1)
        shape = np.column_stack([
            sep, exp, sd, ent, mu3, mu4, m1a, med - exp,
            Mask.sum(1), Mp.sum(1),
            CENTRES[o1] - exp, np.where(has2, CENTRES[o2] - exp, 0.0),
            prob[np.arange(len(o1)), o1], np.where(has2, prob[np.arange(len(o1)), o2], 0.0),
        ])
        F.append(np.column_stack([shape, prob]))
        Y.append(d["dtrue"] - exp)                       # ORACLE label: the signed correction
        T.append(np.full(len(exp), ti))
        FO.append(np.full(len(exp), int(d["fold"])))
        # ORACLE nearest mode + matched min-of-N null, same construction as P2b
        rng = SD.stable_rng(p, "D2")
        pos = np.where(Mask, CENTRES[None], np.nan)
        dd = np.abs(pos - d["dtrue"][:, None])
        hasm = Mask.any(1)
        best = np.where(hasm, CENTRES[np.nanargmin(np.where(np.isnan(dd), np.inf, dd), 1)], exp)
        off = np.abs(pos - exp[:, None])
        cand = exp[:, None] + rng.choice([-1.0, 1.0], size=off.shape) * off
        cd = np.where(np.isnan(cand), np.inf, np.abs(cand - d["dtrue"][:, None]))
        nb = np.where(hasm, cand[np.arange(len(exp)), np.argmin(cd, 1)], exp)
        EXTRA.append(np.column_stack([
            exp, d["dtrue"], deb[int(d["fold"])](exp, sep), best, nb,
            Mask.sum(1), sep, sd]))
    return (np.concatenate(F), np.concatenate(Y), np.concatenate(T),
            np.concatenate(FO), np.concatenate(EXTRA), pdbs, folds)


def ridge(X, y, lam=10.0):
    Xb = np.column_stack([X, np.ones(len(X))])
    A = Xb.T @ Xb + lam * np.eye(Xb.shape[1])
    A[-1, -1] -= lam
    return np.linalg.solve(A, Xb.T @ y)


def apply_ridge(w, X):
    return np.column_stack([X, np.ones(len(X))]) @ w


def main():
    X, y, tid, fo, E, pdbs, folds = build()
    exp, dtrue, dep, best, nullb, nmode, sep, sd = [E[:, k] for k in range(8)]
    rng = SD.stable_rng("D", "P6")
    ntar = len(pdbs)
    mu, sg = X.mean(0), X.std(0) + 1e-9
    Xn = (X - mu) / sg
    noshape = [0, 1, 2]                              # sep, mean, sd only

    def lfo(cols, model="ridge", perm=False, oracle=False):
        pred = np.zeros(len(y))
        for f in sorted(set(fo)):
            te = fo == f
            tr = ~te if not oracle else te
            yt = y[tr].copy()
            if perm:
                yt = rng.permutation(yt)
            if model == "ridge":
                w = ridge(Xn[tr][:, cols], yt)
                pred[te] = apply_ridge(w, Xn[te][:, cols])
            else:
                from sklearn.ensemble import HistGradientBoostingRegressor as G
                g = G(max_iter=200, max_depth=4, learning_rate=0.06, random_state=0)
                g.fit(Xn[tr][:, cols], yt)
                pred[te] = g.predict(Xn[te][:, cols])
        return exp + pred

    def loto(cols):
        """Leave-one-TARGET-out ridge.  Separates 'the relationship does not transfer across
        folds' from 'it does not transfer at all'."""
        pred = np.zeros(len(y))
        for t in range(ntar):
            te = tid == t
            w = ridge(Xn[~te][:, cols], y[~te])
            pred[te] = apply_ridge(w, Xn[te][:, cols])
        return exp + pred

    def within_fold_split(cols):
        """Fit on half the targets of the HELD-OUT fold, predict the other half.  A genuine
        out-of-sample ceiling matched to the fold's own distribution -- unlike the in-fold
        arms below, which are IN-SAMPLE fits and bound the class only loosely."""
        from sklearn.ensemble import HistGradientBoostingRegressor as G
        pred = np.zeros(len(y))
        for f in sorted(set(fo)):
            ts = sorted(set(tid[fo == f]))
            h = set(ts[::2])
            a = np.isin(tid, list(h)) & (fo == f)
            b = (fo == f) & ~a
            for tr, te in ((a, b), (b, a)):
                if tr.sum() < 50 or te.sum() == 0:
                    continue
                g = G(max_iter=200, max_depth=4, learning_rate=0.06, random_state=0)
                g.fit(Xn[tr][:, cols], y[tr])
                pred[te] = g.predict(Xn[te][:, cols])
        return exp + pred

    allc = list(range(X.shape[1]))
    aims = {
        "mean_raw": exp,
        "deployed (shipped debias)": dep,
        "lfo_noshape ridge": lfo(noshape),
        "lfo_shape ridge": lfo(allc),
        "lfo_shape GBM": lfo(allc, "gbm"),
        "lfo_perm ridge (ZERO-INFO)": lfo(allc, perm=True),
        "loto_shape ridge": loto(allc),
        "withinfold-split GBM": within_fold_split(allc),
        "ORACLE_infold ridge (IN-SAMPLE)": lfo(allc, oracle=True),
        "ORACLE_infold GBM (IN-SAMPLE)": lfo(allc, "gbm", oracle=True),
        "ORACLE_bestmode": best,
        "minofN_null": nullb,
    }

    def per_target(v, mask=None):
        m = np.ones(len(v), bool) if mask is None else mask
        return np.asarray([np.abs(v - dtrue)[m & (tid == t)].mean()
                           if (m & (tid == t)).any() else np.nan for t in range(ntar)])

    def boot(dif, B=4000):
        dif = dif[np.isfinite(dif)]
        k = len(dif)
        b = dif[rng.integers(0, k, size=(B, k))].mean(1)
        return float(dif.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))

    for lab, mask in [("ALL PAIRS", None), ("MULTIMODAL PAIRS ONLY", nmode >= 2)]:
        base = per_target(dep, mask)
        print(f"\nP6 -- {lab}:  mean |aim - true| per target, vs the DEPLOYED debiased mean")
        print(f"  {'arm':<30}{'MAE':>8}{'vs deployed [95% CI]':>28}{'W/L':>10}")
        for k, v in aims.items():
            pt = per_target(v, mask)
            m_, lo, hi = boot(pt - base)
            ok = np.isfinite(pt - base)
            print(f"  {k:<30}{np.nanmean(pt):>8.3f}   {m_:+.3f} [{lo:+.3f},{hi:+.3f}]"
                  f"{int(((pt - base) < 0)[ok].sum()):>6}/{int(((pt - base) > 0)[ok].sum())}")
        d1 = per_target(aims["ORACLE_bestmode"], mask) - per_target(aims["minofN_null"], mask)
        m_, lo, hi = boot(d1)
        print(f"  mode information (ORACLE - min-of-N null)   {m_:+.3f} [{lo:+.3f},{hi:+.3f}]")

    # ---------------------------------------------------------------- P2a robustness
    print("\nP2a ROBUSTNESS -- non-parametric nearest-neighbour match on (sep exact, sd nearest),")
    print("within target, multimodal pair -> unimodal pair.  |mean - true| difference.")
    multi = nmode >= 2
    e = np.abs(exp - dtrue)
    per_t = []
    for t in range(ntar):
        m = tid == t
        idx = np.flatnonzero(m)
        A = idx[multi[idx]]
        B = idx[~multi[idx]]
        if len(A) == 0 or len(B) == 0:
            per_t.append(np.nan)
            continue
        ds = []
        for a in A:
            cand = B[sep[B] == sep[a]]
            if not len(cand):
                continue
            b = cand[np.argmin(np.abs(sd[cand] - sd[a]))]
            if abs(sd[b] - sd[a]) > 0.25:               # refuse a bad match
                continue
            ds.append(e[a] - e[b])
        per_t.append(float(np.mean(ds)) if ds else np.nan)
    per_t = np.asarray(per_t)
    m_, lo, hi = boot(per_t)
    ok = np.isfinite(per_t)
    print(f"  matched multi - uni = {m_:+.3f} [{lo:+.3f},{hi:+.3f}]   n_t={int(ok.sum())}  "
          f"W/L {int((per_t[ok] < 0).sum())}/{int((per_t[ok] > 0).sum())}")

    json.dump({"n_targets": ntar, "n_pairs": int(len(y))},
              open(os.path.join(OUT, "meta.json"), "w"))
    open(os.path.join(OUT, "COMPLETE"), "w").write(f"n={ntar}\n")


if __name__ == "__main__":
    main()
