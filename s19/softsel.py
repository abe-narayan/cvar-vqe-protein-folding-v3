"""s19/softsel.py -- LEVER (2): DO NOT CONSUME THE ARGMIN.

WHY THIS EXISTS.  It is the only live lever in Sprint 19 with nobody on it, and three separate
results now point at it:

  * L12 (AGENT C, EXACT).  Gate damage decomposes as
        readout^2 = ||b_pool||^2 + 2||b_pool||*ALIGN + ||delta||^2
    and **score gates raise ALIGN while score-free gates do not** (+0.0209 [+0.0004, +0.0402] vs
    -0.0048 [-0.0124, +0.0020]).  A HARD top-m cut is the most score-ordered operation in the
    pipeline.
  * L11 (AGENT A).  The most harmful shared error component is the one aligned with the
    RETRIEVAL POOL, and a gate that prefers pool-typical geometry selects toward it.
  * L6 (AGENT D).  The fit's argmin is exactly where the in-manifold error gets realised.

And the recorded operator law says the terminal reads the set MEAN, not the set best.  **So the
pipeline applies a hard, score-ordered selection and then consumes a mean.**  That is a
mismatch at the joint, and it has never been tested.

    THE QUESTION: if the top-m cut is replaced by a SOFT, temperature-weighted barycentre over
    the whole pool -- no candidate ever selected, none ever discarded -- does the emitted
    structure improve?

THE FAMILY.  One knob.  With shipped Bayes-risk scores `s` over the K = 500 pool,

    w_c  proportional to  exp(-(s_c - min s) / T)

    T -> 0        the hard argmin (one candidate)
    T small       approaches the shipped hard top-m behaviour
    T -> infinity uniform over the whole pool = the pool mean

so the sweep contains the incumbent's *behaviour* and both zero-information endpoints as limits
of ONE parameter, which is what makes it interpretable.  Averaging is done in the shipped
operator's own frame: superpose on the medoid, then take the WEIGHTED mean.

ARMS.  All native-free.  The native is read only to score.

    top75         the shipped hard cut, coordinate-averaged        THE INCUMBENT, ~3.048 A
    soft{T}       weighted barycentre over all 500, temperature T  the hypothesis
    poolmean      uniform over all 500                             ZERO-INFORMATION endpoint
    rand75        a random 75 of the pool, averaged                MATCHED-RANDOM to the cut
    hard{m}       hard top-m for m in (25, 150, 300)               the cut-size axis, for contrast

CONTROLS AND WHY EACH.  `poolmean` is the T -> infinity endpoint and the zero-information
reference: if soft weighting never beats it, the scores contribute nothing through this
operator.  `rand75` is the matched-random control for the CUT ITSELF at the shipped size -- the
Sprint-18/19 rule is that a physics or score gate must beat discarding the same NUMBER of
candidates at random.  `hard{m}` separates "the cut is too small" from "the cut is too hard",
which are different claims and have been confused before.

PRE-REGISTRATION, written before the run.

  HYPOTHESIS.  Some finite T beats the shipped hard cut at n = 126, because a soft weighting
  raises ALIGN less than a hard score-ordered cut does.

  PRIMARY OUTCOME.  min over the pre-declared T grid of (soft{T} - top75), paired over targets,
  fold-clustered bootstrap CI.  The T grid is fixed here and is NOT extended after seeing
  results.

  EXPECTED, honestly.  I expect this to FAIL, and I am running it because it is cheap and
  unexamined rather than because I believe it.  L12 says score ORDERING is what damages the
  readout, and a soft weighting is still ordered by the same score -- it is a gentler version of
  the same operation, not a different one.  The strongest version of the lever is "never form a
  point estimate anywhere", which this is not.  A plausible outcome is that the whole T curve is
  flat between the two endpoints and the incumbent's m = 75 is already near-optimal, which would
  be a clean negative.

  FALSIFIER.  No T beats `top75` with a CI excluding zero, OR the best T fails to beat
  `poolmean` (which would mean the scores contribute nothing and any gain is the softening
  alone).

  THE TRAP I AM WATCHING FOR.  A win at large T is a win for AVERAGING MORE CANDIDATES, not for
  soft consumption -- `hard150`/`hard300` price that axis separately, and any claimed soft-vs-hard
  effect must survive comparison to a hard cut of the same EFFECTIVE size.  Effective sample
  size (1 / sum w^2) is recorded per target for exactly this reason.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I             # noqa: E402
from s15 import seed as SD                  # noqa: E402

#: pre-declared and not extended after seeing results
TEMPS = (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0)
HARD_M = (25, 150, 300)


def _wavg(W, w):
    """The shipped operator with weights: superpose on the medoid, then WEIGHTED mean.

    `I.coordinate_average` is the unweighted case of exactly this, and is used unchanged for the
    hard arms so the incumbent is reproduced rather than re-implemented.
    """
    W = np.asarray(W, float)
    P = I.pairwise_rmsd(W)
    b = I.medoid(P)
    S = I.superpose_batch(W, W[b])
    w = np.asarray(w, float)
    w = w / max(w.sum(), 1e-300)
    return (S * w[:, None, None]).sum(0)


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows, t0 = [], time.time()

    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        p = I.pool_idx(u)
        Wpool = np.asarray(u["W"], float)[p]
        rec = I.shipped_record(pdb)
        sub = np.asarray(rec["sub"], int)
        rng = SD.stable_rng(pdb, "s19softsel")

        #: the shipped Bayes-risk score over the WHOLE pool -- native-free, the same object the
        #: production filter ranks by.
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        D = I.pair_dists(Wpool, i, j)
        s = np.asarray(I.shipped_score(dg, D), float)
        order = np.argsort(s)

        e = {"n_pool": int(len(p))}

        #: INCUMBENT -- the shipped hard cut, averaged by the unmodified operator
        avg75, _b = I.coordinate_average(Wpool[sub])
        e["top75"] = float(I.ca_rmsd(avg75, nat))

        #: SOFT -- temperature-weighted barycentre over the entire pool, nothing discarded
        z = s - s.min()
        for T in TEMPS:
            w = np.exp(-z / T)
            key = f"soft{T:g}"
            e[key] = float(I.ca_rmsd(_wavg(Wpool, w), nat))
            wn = w / w.sum()
            e["ess_" + key] = float(1.0 / (wn * wn).sum())

        #: ZERO-INFORMATION endpoint (T -> infinity)
        e["poolmean"] = float(I.ca_rmsd(_wavg(Wpool, np.ones(len(Wpool))), nat))

        #: MATCHED-RANDOM to the cut: discard the same NUMBER at random
        r75 = rng.choice(len(Wpool), size=min(75, len(Wpool)), replace=False)
        avgr, _b = I.coordinate_average(Wpool[r75])
        e["rand75"] = float(I.ca_rmsd(avgr, nat))

        #: the CUT-SIZE axis, so "too small" and "too hard" are not confused
        for m in HARD_M:
            mm = min(m, len(Wpool))
            avgm, _b = I.coordinate_average(Wpool[order[:mm]])
            e[f"hard{m}"] = float(I.ca_rmsd(avgm, nat))

        rows.append({"pdb": pdb, "n": n, "fold": fold, **e})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": False},
                      open(os.path.join(RESULTS, "softsel.json"), "w"))
    json.dump({"rows": rows, "complete": len(rows) == len(tg)},
              open(os.path.join(RESULTS, "softsel.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "softsel.json")))["rows"]
    rng = SD.stable_rng("softsel", "report")
    g = lambda k: np.array([r[k] for r in rows])           # noqa: E731
    base = g("top75")
    print(f"\nn = {len(rows)}.  Baseline = the shipped hard top-75 cut, coordinate-averaged "
          f"({base.mean():.3f} A).  All arms NATIVE-FREE.\n")
    print(f"  {'arm':<14}{'RMSD':>8}{'median':>9}{'vs top75 [95% CI]':>27}{'W/L':>9}{'ESS':>8}")
    print(f"  {'top75 (incumb)':<14}{base.mean():>8.3f}{np.median(base):>9.3f}{'--':>27}"
          f"{'--':>9}{75.0:>8.1f}")
    keys = [f"soft{T:g}" for T in TEMPS] + ["poolmean", "rand75"] + [f"hard{m}" for m in HARD_M]
    for k in keys:
        if k not in rows[0]:
            continue
        v = g(k)
        mu, lo, hi = _boot(v - base, rng)
        w = int((v < base).sum()); l = int((v > base).sum())
        ess = f"{g('ess_' + k).mean():>8.1f}" if ("ess_" + k) in rows[0] else (
            f"{float(len(rows[0]) and rows[0]['n_pool']):>8.1f}" if k == "poolmean" else
            (f"{75.0:>8.1f}" if k == "rand75" else f"{float(k[4:]):>8.1f}"))
        tag = "  CONTROL" if k in ("poolmean", "rand75") else ""
        print(f"  {k:<14}{v.mean():>8.3f}{np.median(v):>9.3f}   "
              f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l}{ess}{tag}")

    best = min((f"soft{T:g}" for T in TEMPS), key=lambda k: g(k).mean())
    mu, lo, hi = _boot(g(best) - base, rng)
    print(f"\n  PRIMARY: best soft arm = {best}   {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")
    mu2, lo2, hi2 = _boot(g(best) - g("poolmean"), rng)
    print(f"    and vs the ZERO-INFORMATION endpoint: {mu2:+.3f} [{lo2:+.3f},{hi2:+.3f}]")
    print("\nREAD.  A win at large T is a win for AVERAGING MORE, not for soft consumption --")
    print("compare it to `hard150`/`hard300` at matched effective size before claiming otherwise.")
    print("If the best soft arm does not beat `poolmean`, the scores contribute nothing through")
    print("this operator and any gain is the softening alone.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
