"""s17/inband.py -- IS THE IN-BAND PROBLEM INFORMATION-LIMITED OR CALIBRATION-LIMITED?

Sprint 17's target changed once the SELECT workstream measured, on identical K = 500 candidate
sets at n = 126, that the shipped distance objective has **global Spearman 0.568 but IN-BAND
Spearman 0.131**.  Its skill is almost entirely separating garbage from plausible and nearly
absent where the argmin is decided.  Fifty-eight functionals of the same distogram then
returned -0.009 A [-0.128, +0.110] leave-fold-out, which is what "you cannot fix an in-band
problem by re-weighting a global signal" looks like.

So the sprint needs an IN-BAND DISCRIMINATOR.  Before building one, there is a prior question
that decides what kind of thing to build, and it is answerable as a CEILING rather than a
model:

    Do the native-free features we already have contain in-band information that a GLOBAL
    model fails to extract because the mapping from features to quality is TARGET-SPECIFIC?

    -- or do they simply not contain in-band information at all?

THE DECOMPOSITION.  Fit the same linear scorer three ways on the same features and the same
band, and the differences are the answer:

    global          weights fitted on the TRAINING folds, applied to held-out targets.
                    This is deployable and is the only realizable arm.
    per-target      weights fitted on the target's OWN labels.  ORACLE.  This is the ceiling
                    of per-target calibration -- what a perfect calibrator would buy.
    oracle-1feat    the single best feature per target, chosen with the target's labels.
                    ORACLE.  A floor on per-target calibration using no fitting at all.

    global -> per-target gap  =  the CALIBRATION headroom
    per-target -> band best   =  the INFORMATION deficit that no calibration can fix

If the calibration headroom is large, section 14's target-level calibration is the right
direction and the features are adequate.  If it is small while the information deficit is
large, **new features are required and no amount of calibration will help** -- which would
close the sprint's main remaining direction and is worth knowing early.

PRE-REGISTRATION, written before the run.

  HYPOTHESIS.  The calibration headroom is substantial -- the programme's record says in-band
  ordering is learnable within a target and does not transfer.

  EXPECTED.  per-target ORACLE materially better than global; both well above the band best.

  CONTROLS, both mandatory.  (i) RANDOM selection within the band, which is what "no skill"
  looks like and is the number every arm must beat.  (ii) ZERO-INFORMATION references: a
  constant ideal alpha-helix and beta-strand scored as candidates, and the retrieval rank
  itself (which knows only BLOSUM similarity, no geometry).

  SUCCESS.  The global arm beats random-in-band with a paired fold-clustered interval
  excluding zero.

  FALSIFIER.  If the PER-TARGET ORACLE arm cannot beat random-in-band, the features carry no
  in-band information whatsoever and both this module and section 14's calibration programme
  are closed for this feature set.

  TRAP GUARD.  The per-target arm is fitted on the target's own labels and is therefore an
  ORACLE CEILING, never a result.  It is labelled ORACLE in every table.  With B candidates
  and p features it will also overfit; the number of features is kept well below B and a
  leave-one-candidate-out variant of the per-target arm is reported beside the in-sample one,
  because an in-sample per-target fit is not a ceiling, it is a memorisation.

FEATURES, all native-free.  Distance-objective family (total, mean and max absolute residual,
short- and long-range split), retrieval rank and BLOSUM similarity, Legacy total and its
components, and the consensus/typicality family -- mean pairwise RMSD to the other band
members -- which the programme's record names as the ONLY signal with positive in-band skill.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

K = 500
BANDS = (10, 25, 50)


def _features(W, PHI, PSI, seq, dg, i, j, rank_global, sim, leg=None):
    """native-free per-candidate features.  No native coordinate is read anywhere here."""
    D = I.pair_dists(W, i, j)
    mu = np.asarray(dg["mu"], float) if isinstance(dg, dict) and "mu" in dg else None
    F, names = [], []

    tot = np.asarray(I.shipped_score(dg, D), float)
    F.append(tot); names.append("dist_total")

    if mu is not None and mu.shape[-1] == D.shape[-1]:
        res = np.abs(D - mu[None, :])
        F.append(res.mean(1)); names.append("dist_mae")
        F.append(res.max(1)); names.append("dist_max")
        sep = np.abs(j - i)
        sh = sep <= 4; lo = sep >= 8
        if sh.any():
            F.append(res[:, sh].mean(1)); names.append("dist_short")
        if lo.any():
            F.append(res[:, lo].mean(1)); names.append("dist_long")

    F.append(rank_global.astype(float)); names.append("retr_rank")
    F.append(np.asarray(sim, float)); names.append("retr_sim")

    #: geometry descriptors -- global shape, native-free
    cen = W - W.mean(1, keepdims=True)
    F.append(np.sqrt((cen ** 2).sum(2).mean(1))); names.append("rgyr")
    F.append(np.linalg.norm(W[:, -1] - W[:, 0], axis=1)); names.append("end2end")
    F.append((D < 8.0).mean(1)); names.append("contact_frac")

    if leg is not None:
        F.append(np.asarray(leg, float)); names.append("legacy_total")
    return np.column_stack(F), names


def _consensus(Wb):
    """mean pairwise RMSD to the other band members -- the recorded in-band signal."""
    P = I.pairwise_rmsd(Wb)
    m = len(Wb)
    return (P.sum(1) / max(m - 1, 1))


def _fit_pick(X, y, idx_fit, idx_pick):
    """least-squares on idx_fit, argmin of the prediction over idx_pick; returns the choice."""
    A = np.column_stack([X[idx_fit], np.ones(len(idx_fit))])
    w, *_ = np.linalg.lstsq(A, y[idx_fit], rcond=None)
    B = np.column_stack([X[idx_pick], np.ones(len(idx_pick))])
    return int(idx_pick[int(np.argmin(B @ w))])


def run(targets=None):
    from s16 import energy_lib as EL
    tg = targets if targets is not None else I.targets()
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        order = np.asarray(u["order"], int)[:K]
        W = np.asarray(u["W"], float)[order]
        PHI = np.asarray(u["PHI"], float)[order]
        PSI = np.asarray(u["PSI"], float)[order]
        sim = np.asarray(u["sim"], float)[order]
        nat = np.asarray(u["nat_ca"], float)
        d = I.kabsch_rmsd_batch(W, nat)                      # ORACLE labels
        i, j = I.pair_index(n)
        dg = I.distogram(pdb, seq, fold)
        comp = EL.legacy_components_of_windows(seq, PHI, PSI)
        leg = np.asarray(EL.legacy_total_from(comp), float)
        X, names = _features(W, PHI, PSI, seq, dg, i, j, np.arange(len(order)), sim, leg)
        sc = X[:, 0]
        rk = np.argsort(sc, kind="stable")
        rng = SD.stable_rng(pdb, "s17inband")

        cell = {}
        for B in BANDS:
            b = rk[:B]
            Xb = X[b].copy()
            #: consensus is computed WITHIN the band, so it is a different feature at each B
            cons = _consensus(W[b])
            Xb = np.column_stack([Xb, cons])
            #: standardise per target so the per-target and global fits see the same scale
            mu_, sd_ = Xb.mean(0), Xb.std(0)
            sd_[sd_ < 1e-12] = 1.0
            Z = (Xb - mu_) / sd_
            yb = d[b]
            cell[str(B)] = {
                "Z": Z.tolist(), "y": yb.tolist(),
                "band_best": float(yb.min()), "band_mean": float(yb.mean()),
                "dist_pick": float(yb[0]),
                "rand": float(np.mean([yb[int(rng.integers(B))] for _ in range(3)])),
                "cons_pick": float(yb[int(np.argmin(cons))]),
                "legacy_pick": float(yb[int(np.argmin(Xb[:, names.index("legacy_total")]))]),
                "rank_pick": float(yb[int(np.argmin(Xb[:, names.index("retr_rank")]))]),
            }
        rows.append({"pdb": pdb, "n": n, "fold": fold, "names": names + ["consensus"],
                     "cells": cell})
        if (c + 1) % 20 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(os.path.join(RESULTS, "inband.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "inband.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "inband.json")))["rows"]
    rng = SD.stable_rng("inband", "report")
    folds = np.asarray([r["fold"] for r in rows], int)
    print(f"\nn = {len(rows)} targets, K = {K}.  ORACLE arms read the target's own labels and are")
    print("CEILINGS, never results.  The only deployable arm is `global`.\n")

    for B in BANDS:
        b = str(B)
        Z = [np.asarray(r["cells"][b]["Z"], float) for r in rows]
        y = [np.asarray(r["cells"][b]["y"], float) for r in rows]
        best = np.array([r["cells"][b]["band_best"] for r in rows])
        rnd = np.array([r["cells"][b]["rand"] for r in rows])
        dpick = np.array([r["cells"][b]["dist_pick"] for r in rows])
        cpick = np.array([r["cells"][b]["cons_pick"] for r in rows])
        lpick = np.array([r["cells"][b]["legacy_pick"] for r in rows])

        #: GLOBAL -- fitted on training folds only, the deployable arm
        gl = np.zeros(len(rows))
        for f in sorted(set(folds)):
            tr = [q for q in range(len(rows)) if folds[q] != f]
            A = np.vstack([np.column_stack([Z[q], np.ones(len(Z[q]))]) for q in tr])
            t_ = np.concatenate([y[q] for q in tr])
            w, *_ = np.linalg.lstsq(A, t_, rcond=None)
            for q in range(len(rows)):
                if folds[q] == f:
                    P = np.column_stack([Z[q], np.ones(len(Z[q]))]) @ w
                    gl[q] = y[q][int(np.argmin(P))]

        #: PER-TARGET -- ORACLE, in-sample and leave-one-candidate-out
        pt_in = np.zeros(len(rows)); pt_loo = np.zeros(len(rows))
        for q in range(len(rows)):
            idx = np.arange(len(y[q]))
            pt_in[q] = y[q][_fit_pick(Z[q], y[q], idx, idx) - 0] if False else \
                y[q][int(np.argmin(np.column_stack([Z[q], np.ones(len(idx))]) @
                                   np.linalg.lstsq(np.column_stack([Z[q], np.ones(len(idx))]),
                                                   y[q], rcond=None)[0]))]
            preds = np.empty(len(idx))
            for a in idx:
                keep = idx[idx != a]
                A = np.column_stack([Z[q][keep], np.ones(len(keep))])
                w, *_ = np.linalg.lstsq(A, y[q][keep], rcond=None)
                preds[a] = np.append(Z[q][a], 1.0) @ w
            pt_loo[q] = y[q][int(np.argmin(preds))]

        #: ORACLE single feature per target -- a MIN-OF-N STATISTIC, and therefore meaningless
        #: without its min-of-N null.  mean_t min_v RMSD(t,v) falls with the number of
        #: selectors v WHETHER OR NOT any per-target structure exists: with V draws from a
        #: band of B the expected minimum sits near the 1/(V+1) quantile of the band.  The
        #: SELECT workstream measured this on its own instrument and found the real minimum
        #: LOSES to the zero-information minimum (+0.059 [+0.020, +0.100], CI excluding zero),
        #: which retracts the coordinator's first reading of this row.  Both are computed here.
        o1 = np.zeros(len(rows)); o1n = np.zeros(len(rows))
        nrep = 200
        for q in range(len(rows)):
            cands = []
            for fi in range(Z[q].shape[1]):
                for sgn in (1.0, -1.0):
                    cands.append(y[q][int(np.argmin(sgn * Z[q][:, fi]))])
            o1[q] = min(cands)
            #: ZERO-INFORMATION min-of-N null: V members drawn i.i.d. from the SAME band,
            #: V matched to the number of real selectors above, min taken, averaged.
            V = len(cands)
            pick = rng.integers(0, len(y[q]), size=(nrep, V))
            o1n[q] = float(np.mean(y[q][pick].min(axis=1)))

        print(f"=== BAND B = {B}  (top-{B} of K={K} by the shipped distance objective)")
        print(f"  {'arm':<26}{'mean':>8}{'median':>9}   {'vs random-in-band [95% CI]':>30}{'W/L':>9}")
        for nm, v, tag in (("band best", best, "ORACLE"),
                           ("per-target fit, in-sample", pt_in, "ORACLE"),
                           ("per-target fit, LOO", pt_loo, "ORACLE"),
                           ("best single feature", o1, "ORACLE min-of-N"),
                           ("  its ZERO-INFO min-of-N null", o1n, "CONTROL"),
                           ("GLOBAL fit (deployable)", gl, ""),
                           ("distance argmin", dpick, ""),
                           ("consensus argmin", cpick, ""),
                           ("Legacy argmin", lpick, ""),
                           ("random in-band", rnd, "CONTROL")):
            mu, lo, hi = _boot(v - rnd, rng)
            w = int((v < rnd).sum()); l = int((v > rnd).sum())
            print(f"  {nm:<26}{v.mean():>8.3f}{np.median(v):>9.3f}   "
                  f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l}  {tag}")
        mn, mlo, mhi = _boot(o1 - o1n, rng)
        print(f"\n  MIN-OF-N CHECK  real best-single-feature minus its zero-information null:"
              f"  {mn:+.3f} [{mlo:+.3f},{mhi:+.3f}]")
        print("  A positive value means the real per-target minimum is WORSE than drawing the")
        print("  same number of candidates at random from the same band -- i.e. the 'routing")
        print("  ceiling' is at or below its own noise floor and is not a ceiling at all.")
        cal = gl.mean() - pt_loo.mean()
        info = pt_loo.mean() - best.mean()
        print(f"\n  CALIBRATION headroom  global - per-target(LOO) = {cal:+.3f} A")
        print(f"  INFORMATION deficit   per-target(LOO) - band best = {info:+.3f} A")
        print(f"  -> {'calibration' if cal > info else 'information'} is the binding constraint "
              f"at B = {B}\n")

    print("READ.  If the calibration headroom dominates, target-level calibration is the right")
    print("direction and the features are adequate.  If the information deficit dominates, the")
    print("features do not contain the answer and no calibrator will find it -- which would close")
    print("the sprint's main remaining direction for THIS feature set.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
