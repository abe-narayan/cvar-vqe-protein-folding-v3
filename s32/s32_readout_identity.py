#!/usr/bin/env python
"""s32/s32_readout_identity.py -- the readout identity and its convexity, with a control.

For any weights with sum(w) = 1, candidates W_x and target t:

    || sum_x w_x W_x - t ||^2  =  <w, a>  -  0.5 * w' B w

    a_x  = ||W_x - t||^2     ORACLE      (per-candidate quality; the only unknown)
    B_xy = ||W_x - W_y||^2   native-free (free, computable at inference)

B is a Euclidean squared-distance matrix, hence conditionally negative semidefinite on
{sum v = 0}: for such v, v'Bv = -2||sum_x v_x u_x||^2 <= 0 with u_x = W_x - t. So w'Bw is
CONCAVE on the simplex, -0.5 w'Bw is CONVEX there, and the whole program is linear + convex
= convex, with no hyperparameter.

Reading the two terms: minimising wants LOW <w,a> (individually good candidates) and HIGH
w'Bw (maximally spread candidates -- the variance-cancellation term). **Spread is rewarded,
not penalised.** That is why quality-blind dispersion maximisation, which is this program
with a = const, picks garbage: S31 measured it at +0.1436 A, 1.22x MDE, WORSE.

WHY THIS FILE EXISTS AS A TEST RATHER THAN A COMMENT
The coordinator's first version of the convexity check was VACUOUS: it initialised the
accumulator at the pass threshold (0.0) and took a max over quantities that are always
negative, so it could not fail. That is contract rule 5 -- "a verification must be able to
fail" -- violated by the person who wrote the rule, one hour after writing it. This version
reports BOTH ends of the range and runs a POSITIVE CONTROL (a symmetric non-negative matrix
that is not a squared-distance matrix) which the test must catch.

    python s32/s32_readout_identity.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "s32_readout_identity.json")
SEED = 320001


def sdm(W):
    n = W.shape[0]
    D = W.reshape(n, -1)
    return ((D[:, None, :] - D[None, :, :]) ** 2).sum(-1)


def cnsd(B, rng, trials=600):
    """(max, min) of normalised v'Bv over {sum v = 0}. <= 0 iff conditionally neg. semidef."""
    n = B.shape[0]
    scale = np.abs(B).max() or 1.0
    vals = []
    for _ in range(trials):
        v = rng.normal(size=n)
        v -= v.mean()
        vals.append(float((v @ B @ v) / (scale * (v @ v))))
    return max(vals), min(vals)


def id_err(W, t, rng, trials=400):
    n = W.shape[0]
    a = ((W - t) ** 2).sum(axis=(1, 2))
    B = sdm(W)
    worst = 0.0
    for _ in range(trials):
        w = rng.dirichlet(np.ones(n) * rng.uniform(0.15, 3.0))
        lhs = float((((w[:, None, None] * W).sum(0) - t) ** 2).sum())
        worst = max(worst, abs(lhs - (w @ a - 0.5 * w @ B @ w)) / max(1.0, abs(lhs)))
    return float(worst)


def main():
    rng = np.random.default_rng(SEED)
    n, L = 12, 14
    base = rng.normal(size=(L, 3)) * 3.0
    t = rng.normal(size=(L, 3)) * 3.0
    W = base[None] + rng.normal(size=(n, L, 3)) * 1.2
    Wd = W.copy()
    Wd[3] = Wd[0]; Wd[7] = Wd[0]; Wd[9] = Wd[1]          # duplicates, as real pools have
    Wc = np.stack([base + k * 0.4 * rng.normal(size=(L, 3)) for k in np.linspace(0, 1, n)])
    Ws = np.repeat(W[:1], n, axis=0)                      # every candidate identical: B = 0

    res = {"seed": SEED, "regimes": {}}
    print("%-34s %11s   %11s %11s" % ("regime", "id err", "max v.Bv", "min v.Bv"))
    for nm, X in (("generic", W), ("duplicates", Wd), ("near_collinear", Wc), ("all_identical", Ws)):
        e = id_err(X, t, rng)
        hi, lo = cnsd(sdm(X), rng)
        res["regimes"][nm] = {"identity_rel_err": e, "max_vBv": hi, "min_vBv": lo,
                              "identity_ok": bool(e < 1e-12), "cnsd_ok": bool(hi <= 1e-12)}
        print("%-34s %11.2e   %+11.3e %+11.3e" % (nm, e, hi, lo))

    Bf = np.abs(rng.normal(size=(n, n)))
    Bf = Bf + Bf.T
    np.fill_diagonal(Bf, 0.0)
    hi, lo = cnsd(Bf, rng)
    res["positive_control"] = {"max_vBv": hi, "min_vBv": lo, "caught": bool(hi > 1e-12)}
    print("%-34s %11s   %+11.3e %+11.3e" % ("POSITIVE CONTROL (not a SDM)", "-", hi, lo))

    res["identity_holds_all_regimes"] = all(v["identity_ok"] for v in res["regimes"].values())
    res["cnsd_holds_all_regimes"] = all(v["cnsd_ok"] for v in res["regimes"].values())
    res["selftest_caught_control"] = res["positive_control"]["caught"]
    res["CONCLUSION"] = ("<w,a> - 0.5 w'Bw is linear + convex on the simplex; the program is "
                         "convex and tuning-free. Spread is REWARDED, not penalised.")
    print()
    for k in ("identity_holds_all_regimes", "cnsd_holds_all_regimes", "selftest_caught_control"):
        print("%-42s %s" % (k, res[k]))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(res, fh, indent=1)
    os.replace(tmp, OUT)
    print("\nwrote", OUT)
    return 0 if (res["identity_holds_all_regimes"] and res["cnsd_holds_all_regimes"]
                 and res["selftest_caught_control"]) else 1


if __name__ == "__main__":
    sys.exit(main())
