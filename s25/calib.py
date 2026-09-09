"""s25/calib.py -- IS THE DISTANCE POSTERIOR CALIBRATED?  THE DIAGNOSTIC THAT DECIDES PHASE I.

WHY THIS RUNS FIRST.  Sprint 24's prior ladder (s24/priorladder.py, n=126) measured the endpoint's
sensitivity to prior quality as a transfer function: moving the posterior's probability mass toward
the native by gamma gives -2.1496 A per unit gamma at the origin, concave, and **gamma = 0.0225
reaches 3.0000 A**.  That is an ORACLE result -- it prices what a better prior would BUY and says
nothing about obtainability.

The endgame directive's SS10 asks the only question that can now change the endpoint:
**is the first couple of percent of that ladder REACHABLE from information already available?**

A miscalibrated posterior is the cheapest possible source of such an improvement, because
recalibration needs no new information -- only a correction estimated on OTHER targets.  This file
does not fit anything.  It measures whether there is anything to fit.

WHAT IS MEASURED, per pair, over the score's own pair set (min_sep=2):

    z          = (D_native - expected) / sd            calibrated => mean 0, sd 1
    coverage   = P(D_native inside the posterior's central 50% / 90% interval)
    nll        = -log P(the bin D_native actually falls in)
    modality   = number of local maxima in the 17-bin posterior
    signed     = D_native - expected                   stratified by |i-j| and by chain length

STRATIFICATION IS THE POINT.  A globally calibrated posterior can still be badly miscalibrated
per residue-separation band, and separation is the axis on which a Ca distogram's information
content varies most (the bins themselves are 0.5 A wide at 4 A and 3.5 A wide at 21 A).  If
miscalibration is concentrated in a band, a few-parameter correction fitted cross-target can reach
it; if the posterior is already calibrated everywhere, this whole direction is closed and Phase I
must look elsewhere.

NATIVE USE, STATED PLAINLY.  Native distances are read here for DIAGNOSIS ONLY.  Nothing is fitted,
nothing is selected, and no production object is modified by this file.  Any correction that comes
out of it must be fitted on TRAINING FOLDS ONLY and applied to held-out folds, and must be reported
that way.  Note that the 5 distogram fold models already each saw ~100 of the other 125 dev natives
(s24 Workstream A), so a cross-target correction is strictly LESS exposed than the predictor it
corrects -- but that is a reason to state the discipline, not to relax it.

OPERATOR FORKS, per the standing rule 0.  Diagnostic rather than directional, so the falsifier is
stated in place of a treatment/control pair.  NOTE: enumerated by the coordinator, who has a stake;
the audit lane should re-enumerate.

    functional     DECLARED the shipped posterior `prob` over the shipped `centres`, and `expected`
                   and `sd` exactly as the artefact publishes them.  NOT TAKEN a re-derived mean or
                   sd from `prob` (reported as a consistency check instead, so any disagreement
                   between the published moments and the posterior is visible rather than assumed
                   away).
    basis          DECLARED distance space over the score's own pair set, min_sep=2 -- the space the
                   functional consumes.  NOT TAKEN Cartesian space, which the prior does not have.
    readout        DECLARED per-pair statistics aggregated per target THEN averaged over targets, so
                   long targets do not dominate.  NOT TAKEN pooling all pairs of all targets, which
                   weights by n^2.
    normalisation  DECLARED z by the posterior's OWN published sd.  NOT TAKEN normalising by a fitted
                   or pooled sd, which would hide exactly the miscalibration under test.
    null           DECLARED perfect calibration: z has mean 0 and sd 1, coverage matches nominal.
                   NOT TAKEN a comparison against another predictor, which does not exist here.
    THE LABEL      DECLARED calibration statistics (z, coverage, NLL) AND the signed bias, reported
                   separately, because a posterior can be correctly centred and wrongly scaled or
                   vice versa and the two have different fixes.  NOT TAKEN a single scalar score.

  Falsifier for the whole PHASE-I direction: if z has mean ~0 and sd ~1 uniformly across separation
  bands and lengths, and coverage matches nominal, then the posterior is already calibrated, there is
  no free improvement in it, and recalibration cannot deliver the ladder's first two percent.
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
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402

SEPBANDS = [(2, 2), (3, 3), (4, 5), (6, 8), (9, 15)]


def _save(o, name="calib.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _bin_edges(c):
    """Edges implied by the published bin centres: midpoints, with the outer edges extrapolated."""
    c = np.asarray(c, float)
    mid = 0.5 * (c[1:] + c[:-1])
    lo = c[0] - (mid[0] - c[0]); hi = c[-1] + (c[-1] - mid[-1])
    return np.concatenate([[lo], mid, [hi]])


def run():
    tg = I.targets()
    rows = []
    print("targets: %d" % len(tg), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float); n_res = int(u["n"])
        i, j = I.pair_index(n_res)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        P = np.asarray(dg["prob"], float)                  # (npairs, 17)
        C = np.asarray(dg["centres"], float)
        exp = np.asarray(dg["expected"], float)
        sd = np.maximum(np.asarray(dg["sd"], float), 1e-6)
        Dt = I.pair_dists(nat[None], i, j)[0]
        sep = (j - i).astype(int)

        #: consistency: do the published moments match the posterior they came from?
        m_re = (P * C[None, :]).sum(1)
        v_re = (P * (C[None, :] - m_re[:, None]) ** 2).sum(1)
        sd_re = np.sqrt(np.maximum(v_re, 0.0))

        z = (Dt - exp) / sd
        edges = _bin_edges(C)
        b_true = np.clip(np.digitize(Dt, edges) - 1, 0, len(C) - 1)
        nll = -np.log(np.maximum(P[np.arange(len(P)), b_true], 1e-12))

        #: coverage from the discrete posterior CDF, not from a Gaussian assumption
        cdf = np.cumsum(P, 1)
        q = lambda a: C[np.argmax(cdf >= a, axis=1)]        # noqa: E731
        lo50, hi50, lo90, hi90 = q(0.25), q(0.75), q(0.05), q(0.95)
        cov50 = ((Dt >= lo50) & (Dt <= hi50)).astype(float)
        cov90 = ((Dt >= lo90) & (Dt <= hi90)).astype(float)

        #: multimodality of the 17-bin posterior
        pad = np.pad(P, ((0, 0), (1, 1)), constant_values=0.0)
        peaks = ((pad[:, 1:-1] > pad[:, :-2]) & (pad[:, 1:-1] >= pad[:, 2:]) &
                 (P > 0.02)).sum(1)

        row = {"pdb": pdb, "n": n_res, "fold": int(u["fold"]), "npairs": int(len(Dt)),
               "z_mean": float(z.mean()), "z_sd": float(z.std(ddof=1)) if len(z) > 1 else 0.0,
               "z_absmean": float(np.abs(z).mean()),
               "signed": float((Dt - exp).mean()), "abs_err": float(np.abs(Dt - exp).mean()),
               "sd_mean": float(sd.mean()), "nll": float(nll.mean()),
               "cov50": float(cov50.mean()), "cov90": float(cov90.mean()),
               "multimodal_frac": float((peaks >= 2).mean()),
               "sd_published_vs_posterior": float(np.abs(sd - sd_re).max()),
               "by_sep": {}}
        for a, b in SEPBANDS:
            m = (sep >= a) & (sep <= b)
            if m.sum() == 0:
                continue
            row["by_sep"]["%d-%d" % (a, b)] = {
                "n": int(m.sum()), "z_mean": float(z[m].mean()),
                "z_sd": float(z[m].std(ddof=1)) if m.sum() > 1 else 0.0,
                "signed": float((Dt - exp)[m].mean()), "abs_err": float(np.abs(Dt - exp)[m].mean()),
                "sd_mean": float(sd[m].mean()), "cov90": float(cov90[m].mean()),
                "nll": float(nll[m].mean()),
            }
        rows.append(row)
        if (c_i + 1) % 25 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("z_mean", "z_sd", "cov90", "nll", "by_sep")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "sepbands": SEPBANDS})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "calib.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    print("\nn = %d targets, %d pairs total.  Native distances read for DIAGNOSIS ONLY.\n"
          % (len(rows), int(g("npairs").sum())))
    print("  published sd vs sd recomputed from the posterior: max abs diff %.4f A"
          % g("sd_published_vs_posterior").max())

    print("\n  GLOBAL CALIBRATION  (perfect = z_mean 0.00, z_sd 1.00, cov50 0.50, cov90 0.90)")
    print("    z mean        %+.4f      <- >0 means the posterior predicts distances TOO SHORT"
          % g("z_mean").mean())
    print("    z sd           %.4f      <- >1 means the posterior is OVER-CONFIDENT" % g("z_sd").mean())
    print("    coverage 50%%   %.4f" % g("cov50").mean())
    print("    coverage 90%%   %.4f" % g("cov90").mean())
    print("    signed error  %+.4f A    abs error %.4f A    published sd %.4f A"
          % (g("signed").mean(), g("abs_err").mean(), g("sd_mean").mean()))
    print("    mean NLL of the true bin  %.4f" % g("nll").mean())
    print("    pairs with a MULTIMODAL posterior (>=2 peaks above 0.02): %.3f"
          % g("multimodal_frac").mean())

    print("\n  BY RESIDUE SEPARATION  -- the axis a Ca distogram's information varies on")
    print("    %-8s%8s%10s%9s%11s%11s%10s%9s" %
          ("sep", "pairs", "z_mean", "z_sd", "signed", "abs_err", "sd_pub", "cov90"))
    for a, b in SEPBANDS:
        k = "%d-%d" % (a, b)
        sel = [r["by_sep"][k] for r in rows if k in r["by_sep"]]
        if not sel:
            continue
        f = lambda x: float(np.mean([s[x] for s in sel]))     # noqa: E731
        print("    %-8s%8d%+10.4f%9.4f%+11.4f%11.4f%10.4f%9.4f"
              % (k, int(sum(s["n"] for s in sel)), f("z_mean"), f("z_sd"),
                 f("signed"), f("abs_err"), f("sd_mean"), f("cov90")))

    print("\n  BY CHAIN LENGTH")
    print("    %-8s%8s%10s%9s%11s%9s" % ("n", "targets", "z_mean", "z_sd", "abs_err", "cov90"))
    L = g("n").astype(int)
    for lo, hi in ((9, 10), (11, 12), (13, 14), (15, 16)):
        m = (L >= lo) & (L <= hi)
        if m.sum() == 0:
            continue
        print("    %-8s%8d%+10.4f%9.4f%11.4f%9.4f"
              % ("%d-%d" % (lo, hi), int(m.sum()), g("z_mean")[m].mean(), g("z_sd")[m].mean(),
                 g("abs_err")[m].mean(), g("cov90")[m].mean()))

    print("\n  Falsifier: z_mean ~ 0 and z_sd ~ 1 uniformly, with coverage at nominal, closes the")
    print("  recalibration route -- there would be no free improvement in the posterior to reach.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
