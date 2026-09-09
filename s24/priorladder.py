"""s24/priorladder.py -- B-3.  THE PRIOR-ERROR ATTRIBUTION LADDER.  How much of the incumbent's
3.0483 A is the distogram's own error, and what is the shape of the transfer function?

ORACLE BY CONSTRUCTION.  Every arm interpolates the prior toward the native's own distances.  This
is a TRANSFER FUNCTION, never a system result, and ORACLE is carried on every appearance.

WHAT IS HELD FIXED.  The same shipped K=500 BLOSUM pool, the same shipped Bayes-risk functional,
the same top-75 rule, the same uniform coordinate average in the medoid frame.  ONLY the prior
changes.  So every difference along the ladder is attributable to the prior and to nothing else.

------------------------------------------------------------------------------------------------
THE TWO HAZARDS THE COORDINATOR NAMED, BOTH ADDRESSED BEFORE THE RUN

(1) CAN THE SHIPPED LOSS BE REBUILT EXACTLY?  YES, BIT-EXACTLY, and it is verified in-process by
`_verify_functional()` on every run rather than asserted here.  `core.predict.Distogram` builds

    risk[p, t] = w[p] * sum_c prob[p, c] * |grid[t] - CENTRES[c]|,     w[p] = shell / (sd+0.5)^g

so the functional is a bilinear form in `prob` with a per-pair weight that is itself a function of
`prob` through `sd`.  Rather than re-implement it -- which is exactly the silent approximation the
coordinator warned against -- this module CONSTRUCTS A GENUINE `Distogram` OBJECT from the modified
`prob` and uses its own `_risk`.  Reconstruction of the shipped artefact from the shipped `prob` is
verified at max abs error 0.0.  Note the consequence, which is a feature: `w` is recomputed
self-consistently, so a sharper prior legitimately reweights its own pairs, exactly as a genuinely
better distogram would.  `MASSFIXW` holds `w` at its real value to isolate that channel.

(2) METRIC REALISABILITY.  Interpolating in DISTANCE space, `Dhat(g) = (1-g) Dhat + g Dnat`, can
produce a matrix no structure attains, so the risk table would be scoring an unattainable target.
**This module therefore does not interpolate in distance space at all.**  It interpolates the
PROBABILITY MASS, and every arm's g=1 endpoint is the native's own bin -- realisable by
construction, because the native realises it.  The distance-space arm is reported only as `DIRAC`,
which places the mass on the bin containing the interpolated distance; it is a proper distribution
over a realisable bin at every g, and it is carried as the COLLAPSE control rather than as the
primary.

Stated plainly: the coordinator's literal "substitute Dhat" arm is NOT expressible in the shipped
functional, because that functional consumes a DISTRIBUTION and not a point estimate.  Writing it
as a point-estimate loss would be a different functional and would make the ladder uninterpretable.
Three proper-distribution arms are used instead, and they separate the two things a "better prior"
could mean.

------------------------------------------------------------------------------------------------
THE ARMS -- and why there are three, not one.  "A better distogram" is ambiguous between moving the
prediction's CENTRE and changing its CONFIDENCE.  These separate them:

    MASS(g)   prob <- (1-g) * prob + g * onehot(native bin)
              moves the centre AND sharpens.  PRIMARY: it is the literal "mix the prior toward the
              truth", and its endpoint is a Dirac on the true bin.
    TILT(g)   prob <- prob * exp(-lam |CENTRES - Dnat|), renormalised, with lam solved by bisection
              so that E[d] hits (1-g) E_real + g Dnat.
              moves the centre with MINIMUM KL, preserving the prior's support and most of its
              shape: "a better-CENTRED prior of the same confidence".
    DIRAC(g)  prob <- onehot(bin of (1-g) E_real + g Dnat)
              moves the centre AND collapses completely: "a confidently interpolated prior".
              The COLLAPSE control -- project memory prices posterior collapse at -2.253 A.
    MASSFIXW(g)  MASS with `w` frozen at the real prior's value.  Isolates the reweighting channel.
    EXACT        score = mean_p w_p |D_cand - D_native|, UNBINNED.
              The true perfect-distance-knowledge endpoint.  The 17 bins are coarse (centres 4.0 to
              25.0), so even a Dirac on the native BIN is not perfect distance knowledge, and
              without this arm the ladder's g=1 would be mistaken for it.  DIFFERENT FUNCTIONAL,
              labelled as such, used only as an endpoint reference.

REPRODUCTION GATE.  `EXACT` and `MASS(1)` are checked against the project's standing "perfect
distance knowledge caps the instrument at ~1.95-2.0 A".  If they do not land there the LADDER is
wrong, not the standing number, and that is how it will be reported.

BETA at every g, defined exactly as `s24/referent.py` defines it:
    beta = <Dc - Dt, Dhat - Dt> / |Dhat - Dt|^2
reported against the arm's OWN prior (the one actually in use) and against the REAL prior (so it is
commensurable with L5's +0.5201).  NOTE THE DEGENERACY: as g -> 1 the denominator goes to zero by
construction, so beta is unstable at the top of the ladder and is reported only where the
denominator is non-degenerate.  Flagged rather than smoothed.

RULE 0 -- SIX OPERATOR FORKS, each naming the alternative not taken.

    functional     DECLARED the shipped Bayes-risk score, rebuilt BIT-EXACTLY by constructing a real
                   `core.predict.Distogram` from the modified prob.  NOT TAKEN re-implementing the
                   risk table by hand, or substituting a point-estimate L1 loss, either of which
                   would silently change the functional the ladder claims to measure.
    basis          DECLARED point cloud throughout, native frame for bias quantities.  NOT TAKEN the
                   built-chain basis, which is 0.156 A away and not comparable to the incumbent.
    readout        DECLARED the same top-75 uniform coordinate average in the medoid frame at every
                   g.  NOT TAKEN re-optimising m as the prior improves, which would let the readout
                   absorb the prior's effect and flatten the transfer function.
    normalisation  DECLARED interpolation of PROBABILITY MASS with `w` recomputed self-consistently.
                   NOT TAKEN interpolation in distance space (not realisable, and not expressible in
                   this functional) and NOT TAKEN freezing `w`, which is carried separately as
                   MASSFIXW rather than chosen silently.
    null           DECLARED g=0, which is the incumbent re-derived through this module's own code
                   path and therefore also the module's self-check; plus MASSFIXW and DIRAC as
                   mechanism controls.  NOT TAKEN comparing to the pinned constant alone, which
                   would not catch a defect in this module's own scoring path.
    THE LABEL      DECLARED emitted continuous Ca-RMSD and beta, both at every g, ORACLE throughout.
                   NOT TAKEN reporting the slope alone, which would hide whether the curve is
                   linear, saturating or threshold-like -- the shape is the actual deliverable.

  H  emitted RMSD falls substantially and smoothly with g, so the prior's error is a large and
     incrementally addressable share of the incumbent's 3.0483.
  Falsifier  a flat or threshold-shaped curve at low g.  That would say incremental prior
     improvement buys nothing and only a large jump pays -- which is a DIFFERENT and equally
     decision-relevant answer, and is the one that would redirect the sprint.
"""
from __future__ import annotations

import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s24 import stats_lib as ST           # noqa: E402
import s24.residlib as RL                 # noqa: E402
from core import predict as dgm           # noqa: E402

GAM = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
ARMS = ("MASS", "TILT", "DIRAC", "MASSFIXW")
TOPM = 75
OUT = os.path.join(RL.RES, "priorladder.json")
CEN = np.asarray(dgm.CENTRES, float)


def _verify_functional(pdb="1A13"):
    """Bit-exactness gate, run in-process every time.  If this fails, the ladder is meaningless."""
    u = I.load_univ(pdb)
    dg = I.distogram(pdb)
    d = dgm.Distogram(u["seq"], np.asarray(dg["prob"], float),
                      np.asarray(dg["i"]), np.asarray(dg["j"]))
    err = float(np.abs(np.asarray(dg["risk"], float) - np.asarray(d._risk, float)).max())
    if err != 0.0:
        raise SystemExit("FUNCTIONAL NOT BIT-EXACT (max abs err %g). Stopping rather than "
                         "approximating the loss silently." % err)
    return err


def _tilt(prob, dnat, target):
    """Exponential tilt of `prob` toward `dnat`, lam solved so E[d] == target.  Minimum-KL move."""
    out = np.empty_like(prob)
    for p in range(len(prob)):
        q0 = prob[p]; t = target[p]
        k = np.abs(CEN - dnat[p])
        lo, hi = 0.0, 50.0
        e = lambda lam: float((q0 * np.exp(-lam * k) / (q0 * np.exp(-lam * k)).sum() * CEN).sum())  # noqa: E731
        if abs(e(0.0) - t) < 1e-9 or (e(0.0) - t) * (e(hi) - t) > 0:
            lam = 0.0 if abs(e(0.0) - t) < abs(e(hi) - t) else hi
        else:
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if (e(lo) - t) * (e(mid) - t) <= 0:
                    hi = mid
                else:
                    lo = mid
            lam = 0.5 * (lo + hi)
        w = q0 * np.exp(-lam * k)
        out[p] = w / max(w.sum(), 1e-300)
    return out


def run():
    _verify_functional()
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        n = u["n"]; nat = np.asarray(u["nat_ca"], float); seq = u["seq"]
        idx = I.pool_idx(u); W = np.asarray(u["W"][idx], float)
        dg = I.distogram(pdb)
        pi, pj = np.asarray(dg["i"]), np.asarray(dg["j"])
        prob0 = np.asarray(dg["prob"], float)
        D = I.pair_dists(W, pi, pj)                       # (500, npairs) candidate distances
        Dt = I.pair_dists(nat, pi, pj)                    # native distances  ORACLE
        nb = np.abs(CEN[None, :] - Dt[:, None]).argmin(1)  # native bin per pair  ORACLE
        onehot = np.eye(len(CEN))[nb]
        E0 = (prob0 * CEN[None]).sum(1)
        w_real = dgm.Distogram(seq, prob0, pi, pj).w

        rec = {"pdb": pdb, "n": int(n), "fold": int(u["fold"])}
        rec["P0"] = float(I.ca_rmsd(RL.readout(W[np.argsort(
            I.shipped_score(dg, D), kind="stable")[:TOPM]]), nat))

        def emit_score(sc, tag, Dhat):
            o = np.argsort(sc, kind="stable")[:TOPM]
            C = RL.readout(W[o])
            rec[tag] = float(I.ca_rmsd(C, nat))
            Dc = I.pair_dists(C, pi, pj)
            eC = Dc - Dt
            for nm, dh in (("", Dhat), ("R", E0)):            # own prior, and the REAL prior
                eP = dh - Dt; den = float((eP * eP).sum())
                rec["beta%s_%s" % (nm, tag)] = float((eC * eP).sum() / den) if den > 1e-9 else float("nan")
            rec["den_" + tag] = float(((Dhat - Dt) ** 2).sum())
            rec["ov_" + tag] = float(len(np.intersect1d(o, np.argsort(
                I.shipped_score(dg, D), kind="stable")[:TOPM])) / TOPM)

        for g in GAM:
            tgt = (1 - g) * E0 + g * Dt
            pm = (1 - g) * prob0 + g * onehot
            pt = _tilt(prob0, Dt, tgt) if 0 < g < 1 else (prob0 if g == 0 else onehot)
            pd = np.eye(len(CEN))[np.abs(CEN[None, :] - tgt[:, None]).argmin(1)]
            for arm, pr in (("MASS", pm), ("TILT", pt), ("DIRAC", pd)):
                d = dgm.Distogram(seq, pr, pi, pj)
                emit_score(I.shipped_score({"grid": d.grid, "risk": np.asarray(d._risk, float)}, D),
                           "%s%.1f" % (arm, g), (pr * CEN[None]).sum(1))
            #: MASSFIXW -- same prob as MASS, w frozen at the real prior's value
            dm = dgm.Distogram(seq, pm, pi, pj)
            rk = np.asarray(dm._risk, float) / np.maximum(dm.w[:, None], 1e-12) * w_real[:, None]
            emit_score(I.shipped_score({"grid": dm.grid, "risk": rk}, D),
                       "MASSFIXW%.1f" % g, (pm * CEN[None]).sum(1))

        #: EXACT -- unbinned perfect distance knowledge.  DIFFERENT FUNCTIONAL, endpoint reference.
        emit_score((w_real[None] * np.abs(D - Dt[None])).mean(1), "EXACT", Dt)
        rows.append(rec)
        del u, W, D
        if (c + 1) % 10 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            ST.save_atomic(OUT, {"rows": rows, "gam": GAM, "arms": list(ARMS)},
                           rows=rows, n_expected=len(tg), module_file=__file__)

    need = ["P0", "EXACT"] + ["%s%.1f" % (a, g) for a in ARMS for g in GAM]
    ST.save_atomic(OUT, {"rows": rows, "gam": GAM, "arms": list(ARMS), "topm": TOPM},
                   complete_keys=need, rows=rows, n_expected=len(tg), module_file=__file__)
    report(rows)


def report(rows=None):
    import json
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    g_ = lambda k: np.array([r[k] for r in rows], float)        # noqa: E731
    base = g_("MASS0.0")
    print("\nB-3 / priorladder.  n = %d.  POINT CLOUD.  ORACLE TRANSFER FUNCTION, not a system"
          % len(rows))
    print("result.  Same K=500 pool, same shipped Bayes-risk functional (rebuilt BIT-EXACTLY),")
    print("same top-75, same uniform coordinate average.  ONLY THE PRIOR CHANGES.\n")
    print("  self-check: g=0 through this module's own path = %.4f  (pinned incumbent 3.0483,"
          % base.mean())
    print("              module's own P0 = %.4f)\n" % g_("P0").mean())
    print("  %-6s %9s %9s %9s %9s   %8s %8s" % ("gamma", "MASS", "TILT", "DIRAC", "MASSFIXW",
                                                "betaMASS", "ovMASS"))
    for g in GAM:
        bm = g_("betaR_MASS%.1f" % g)
        print("  %-6.1f %9.4f %9.4f %9.4f %9.4f   %+8.4f %8.3f"
              % (g, g_("MASS%.1f" % g).mean(), g_("TILT%.1f" % g).mean(),
                 g_("DIRAC%.1f" % g).mean(), g_("MASSFIXW%.1f" % g).mean(),
                 np.nanmean(bm), g_("ov_MASS%.1f" % g).mean()))
    print("\n  EXACT (unbinned |D - Dnat|, DIFFERENT FUNCTIONAL, endpoint reference)  %.4f"
          % g_("EXACT").mean())
    print("  betaMASS is measured against the REAL prior at every g, so it is commensurable")
    print("  with L5's +0.5201.  ovMASS = share of the top-75 shared with the incumbent's.")
    print("\n  LOCAL SLOPE  dRMSD/dgamma  (the value of a marginally better prior)")
    for a in ("MASS", "TILT", "DIRAC"):
        v = np.array([g_("%s%.1f" % (a, g)).mean() for g in GAM])
        i1 = min(1, len(v) - 1); i3 = min(3, len(v) - 1)
        print("    %-9s g=0->%.1f %+.4f   0->%.1f %+.4f   full 0->1 %+.4f   curvature %+.4f"
              % (a, GAM[i1], v[i1] - v[0], GAM[i3], v[i3] - v[0], v[-1] - v[0],
                 (v[0] + v[-1]) / 2 - v[len(v) // 2]))
    print("    curvature > 0 = SATURATING (early gamma pays most);  < 0 = THRESHOLD (only a big")
    print("    jump pays).  This shape, not the endpoint, is the deliverable.")
    print("\n  SIGNIFICANCE of each MASS step against g=0")
    for g in GAM[1:]:
        st = RL.stats(g_("MASS%.1f" % g) - base, fold)
        print("    g=%.1f  %+.4f  SE %.4f MDE %.4f  %5.2fx  fold[%+.3f,%+.3f]  %3dW/%3dL  %s"
              % (g, st["mean"], st["se"], st["mde"], st["eff_over_mde"],
                 st["ci_fold"][0], st["ci_fold"][1], st["W"], st["L"], RL.verdict(st)))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
