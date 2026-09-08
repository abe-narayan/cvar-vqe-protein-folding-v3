"""s21/rgsign.py -- THE ONE NATIVE-FREE CHANNEL NEITHER LANE TESTED.

L25 closed the bimodal-containment route on ELEVEN native-free predictors.  Every one of them was a
functional of THE OBJECTIVE'S OWN SCORE DISTRIBUTION, plus window geometry -- and the winner died as
a length artefact.  So what was actually established is narrower than "no native-free signal exists":
it is that **the objective cannot audit itself**.

A different channel exists, is already in the pipeline, and has never been read this way.
`in-band-ordering-is-per-target` records it under a heading that says so explicitly -- *"AND THAT
DIRECTION IS NOT CLOSED -- measured, 126 targets"*:

    the shipped distogram predicts every pair distance, so the RADIUS OF GYRATION follows
    in closed form:   Rg^2 = (1/(2 N^2)) * sum_ij d_ij^2

Length-residualised against the native Rg it scores **0.313 [+0.110,+0.518]**, the retrieval pool's
mean Rg **0.244 [+0.031,+0.444]**, and the incumbent's emitted Rg **0.365 [+0.155,+0.577]** -- all
three CIs excluding zero, against an ORACLE of 0.909.  Real, free, and weak.  That memory's own
verdict is that **the only direction with leverage is anything supplying the PER-TARGET SIGN at
inference: a conditioning signal, not a better objective.**

THE QUESTION THIS FILE ASKS, under WORKSTREAM D's protocol so the answer is comparable to its eleven:

    Does a COMPACTNESS channel predict where the ORACLE-best latent configuration sits in the
    objective's ordering, after length AND target difficulty are partialled out?

WHY A COMPACTNESS CHANNEL MIGHT DO WHAT SCORE FUNCTIONALS CANNOT.  The objective is a distogram
Bayes risk; its score distribution is a summary of how well candidates match the PREDICTED distances.
It has no way to represent the possibility that the PREDICTION ITSELF is wrong on this target.  Rg
computed from the distogram, compared against Rg realised in the pool, is exactly a measure of that
disagreement -- and disagreement between the objective's belief and what the candidate source can
supply is the natural native-free signature of a target where the objective's window is a
concentrated wrong region.

THE SIX PREDICTORS, pre-declared as ONE FAMILY so the best-of-K null is over six and not over one.

    rg_disto        Rg from the distogram's predicted distances        the memory's 0.313 proxy
    rg_pool_mean    mean Rg over the K=500 retrieval pool              the memory's 0.244 proxy
    rg_pool_sd      spread of Rg across the pool                       how much room the source has
    rg_emit         Rg of the incumbent's top-75 coordinate average    the memory's 0.365 proxy
    rg_gap          rg_disto - rg_pool_mean                            THE DISAGREEMENT TERM
    rg_z            (rg_disto - rg_pool_mean) / rg_pool_sd             the same, scaled by the pool

`rg_emit` IS PARTLY CIRCULAR and is labelled so at every appearance -- it is a property of the very
output it would condition, which the source memory flags as a limitation of its 0.365.  It is
carried because excluding it would be choosing which arms to run after knowing what they mean.

OPERATOR FORKS, DECLARED PER BRIEF SS7 RULE 0 -- all six axes, including the one D's failure added.

    functional     DECLARED the shipped Bayes-risk distogram, whose `dhat` IS the prediction being
                   read.  NOT TAKEN the squared functional; it has no `dhat` of its own to read.
    basis          DECLARED the pool's own window coordinates for `rg_pool_*` and `rg_emit`, since
                   Rg is a property of a structure and the shipped pipeline emits windows.  NOT
                   TAKEN the rebuild; L20 prices that difference at +0.016 on this operator.
    readout        DECLARED Spearman (rank), because the claim is ordinal and the distributions are
                   skewed.  NOT TAKEN Pearson, which a single compact outlier would dominate.
    normalisation  DECLARED length-residualised, because chain length drives Rg mechanically
                   (Rg ~ N^0.6 even for a random coil).  NOT TAKEN raw Rg, which would reproduce
                   D's failure exactly.
    null           DECLARED a permutation null on the LABEL with the maximum over all six
                   recomputed inside each permutation (D11).  NOT TAKEN a per-predictor null,
                   which is the wrong null for a family and is the error D11 corrects.
    THE LABEL      DECLARED `ORACLE_rank_pct`, CONTINUOUS.  NOT TAKEN `ORACLE_contained512` --
                   **the binarised label whose threshold is `512/2**n` and therefore
                   LENGTH-DEPENDENT, which manufactured D's 0.776 artefact.**  This axis exists in
                   rule 0 only because that happened, and this file is the first to declare it.

PRE-REGISTRATION.
  Hypothesis    weak or nothing.  The memory's own chain is 0.313 proxy-to-Rg and Rg-to-sign was
                measured on a smaller enumerated set, so composing them ASSUMES linearity and
                independent errors -- the memory says so.  Registered estimate: |rho| < 0.30.
  Falsifier     if the best Rg predictor's partial |rho| exceeds the permuted best-of-6 bar AND its
                own CI excludes zero, a native-free conditioning signal for the bimodal regime
                EXISTS, the 0.4 A is live, and this is the sprint's first positive result.
  Disposition   anything not clearing the permuted bar is NOT DEMONSTRATED, in those words.
  Scope         n = 75 (D's panel, n <= 13), inherited and not to be generalised.
  Promotion     none.  Clearing the bar buys a replication, not a claim.
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

from s12 import instrument as I             # noqa: E402
from s15 import seed as SD                  # noqa: E402

CA_BOND = 3.8046                 #: the trans virtual CA-CA distance; the only |i-j|=1 term
RG_FAMILY = ("rg_disto", "rg_pool_mean", "rg_pool_sd", "rg_emit", "rg_gap", "rg_z")
CHECK = "rg_disto_risk"        #: the undeclared variant, reported but NOT in the family


def _rg2_from_pairs(d2_sum, n):
    """Rg^2 = (1/N^2) * sum_{i<j} d_ij^2.  `d2_sum` is that upper-triangular sum."""
    return d2_sum / float(n * n)


def _rg_of(P):
    """Rg of an (n,3) structure, directly."""
    P = np.asarray(P, float)
    c = P - P.mean(0, keepdims=True)
    return float(np.sqrt((c * c).sum() / len(P)))


def _save(obj):
    path = os.path.join(RESULTS, "rgsign.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def run():
    lr = json.load(open(os.path.join(RESULTS, "d_lrank.json")))["rows"]
    keep = {r["pdb"]: r for r in lr}
    rows = []
    print("targets from d_lrank: %d" % len(keep), flush=True)
    for c, pdb in enumerate(sorted(keep)):
        u = I.load_univ(pdb)
        n = int(u["n"])
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(n)

        #: --- Rg from the DISTOGRAM's own predicted distances.  The distogram supplies |i-j|>=2;
        #: the |i-j|=1 terms are the trans virtual bond, a constant applied identically to every
        #: arm, which is the ideal-geometry approximation the source memory flags.
        #: `dg["expected"]` IS the distogram's posterior-mean predicted distance per pair -- the
        #: certified quantity, not a reconstruction of one.  FORK NOT TAKEN: qb_lib's `dhat`, which
        #: is this debiased by a per-fold separation term; that debias exists to make the OBJECTIVE
        #: well-calibrated, and reading it here would mix the objective's own correction into a
        #: channel whose purpose is to be independent of the objective's scores.
        #: ALSO REPORTED: `rg_disto_risk`, the risk-MINIMISING point estimate, which is what an
        #: undeclared first version of this file used.  Both are carried so the declared choice can
        #: be checked rather than trusted.
        dhat = np.asarray(dg["expected"], float)
        i, j = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
        dhat_risk = np.asarray(dg["grid"], float)[
            np.argmin(np.asarray(dg["risk"], float), axis=1)]
        rg_disto_risk = float(np.sqrt(_rg2_from_pairs(
            float((dhat_risk ** 2).sum() + (n - 1) * CA_BOND ** 2), n)))
        d2 = float((dhat ** 2).sum() + (n - 1) * CA_BOND ** 2)
        rg_disto = float(np.sqrt(_rg2_from_pairs(d2, n)))

        rgs = np.array([_rg_of(w) for w in W], float)
        avg, _b = I.coordinate_average(W[np.argsort(
            np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float), kind="stable")[:75]])
        rows.append({
            "pdb": pdb, "n": n, "fold": int(u["fold"]),
            "rg_disto": rg_disto,
            "rg_disto_risk": rg_disto_risk,
            "rg_pool_mean": float(rgs.mean()),
            "rg_pool_sd": float(rgs.std(ddof=1)),
            "rg_emit": _rg_of(np.asarray(avg, float)),          # PARTLY CIRCULAR
            "rg_gap": float(rg_disto - rgs.mean()),
            "rg_z": float((rg_disto - rgs.mean()) / max(rgs.std(ddof=1), 1e-9)),
            "ORACLE_rank_pct": float(keep[pdb]["ORACLE_rank_pct"]),
            "latent_mean": float(keep[pdb]["latent_mean"]),
        })
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(keep)), flush=True)
    ok = len(rows) == len(keep) and all(
        all(np.isfinite(r[k]) for k in RG_FAMILY) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(keep), "family": list(RG_FAMILY)})
    report(rows)
    return rows


def _rank(x):
    x = np.asarray(x, float)
    o = np.argsort(x, kind="stable")
    r = np.empty(len(x), float)
    r[o] = np.arange(len(x), dtype=float)
    #: average ranks over exact ties, so a tie cannot read the sort order (project rule)
    for v in np.unique(x):
        m = x == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def _partial_spearman(a, b, ctrls):
    """Spearman of a vs b, both residualised on the ranks of every control."""
    ra, rb = _rank(a), _rank(b)
    if ctrls:
        C = np.column_stack([_rank(c) for c in ctrls] + [np.ones(len(ra))])
        ra = ra - C @ np.linalg.lstsq(C, ra, rcond=None)[0]
        rb = rb - C @ np.linalg.lstsq(C, rb, rcond=None)[0]
    else:
        ra = ra - ra.mean(); rb = rb - rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d > 0 else 0.0


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "rgsign.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)       # noqa: E731
    y = g("ORACLE_rank_pct")
    rng = SD.stable_rng("rgsign", "null")
    n_perm = 4000
    print("\nn = %d targets (D's panel).  LABEL: ORACLE_rank_pct, CONTINUOUS." % len(rows))
    print("Family of %d, pre-declared.  Null = permuted label, MAXIMUM over the family.\n"
          % len(RG_FAMILY))

    for lab, ctrls in (("length only", [g("n")]),
                       ("length + DIFFICULTY (binding)", [g("n"), g("latent_mean")])):
        vals = {k: _partial_spearman(g(k), y, ctrls) for k in RG_FAMILY}
        best = max(vals, key=lambda k: abs(vals[k]))
        mx = np.empty(n_perm)
        for t in range(n_perm):
            yp = y[rng.permutation(len(y))]
            mx[t] = max(abs(_partial_spearman(g(k), yp, ctrls)) for k in RG_FAMILY)
        bar = float(np.percentile(mx, 95))
        print("  CONTROL: %s" % lab)
        for k in RG_FAMILY:
            star = "  <- best" if k == best else ""
            circ = "  (PARTLY CIRCULAR)" if k == "rg_emit" else ""
            print("    %-14s partial rho %+.4f%s%s" % (k, vals[k], star, circ))
        #: Bootstrap CI on the BEST arm.  Its point estimate is a MAXIMUM over the family and is
        #: therefore UPWARD-BIASED (D11): the permuted bar handles selection for the SIGNIFICANCE
        #: test, but the magnitude here is not an unbiased estimate of that arm's own value.
        xb = g(best)
        bs = np.array([_partial_spearman(xb[k], y[k], [c[k] for c in ctrls])
                       for k in rng.integers(0, len(y), size=(2000, len(y)))])
        lo, hi = (float(v) for v in np.percentile(bs, [2.5, 97.5]))
        fires = abs(vals[best]) > bar and lo * hi > 0
        print("    permuted best-of-%d bar (95th pct): %.4f   |  best %.4f  CI [%+.3f,%+.3f]  ->  %s"
              % (len(RG_FAMILY), bar, abs(vals[best]), lo, hi,
                 "CLEARS" if fires else "NOT DEMONSTRATED"))
        print("    [check arm, deliberately NOT in the family] %s  partial rho %+.4f"
              % (CHECK, _partial_spearman(g(CHECK), y, ctrls)))
        print("    the best arm is a MAXIMUM OVER 6 and its magnitude is upward-biased.\n")
    print("Registered: |rho| < 0.30 expected; falsifier = best clears the permuted bar with a CI")
    print("excluding zero.  Anything not clearing it is NOT DEMONSTRATED, in those words.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
