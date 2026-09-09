"""s25/temper.py -- THE POSTERIOR IS OVER-CONFIDENT BY 2x. WIDEN IT. ONE GLOBAL PARAMETER, NESTED CV.

WHY.  `s25/calib.py` (n=126, 8,549 pairs) measured the shipped distance posterior against its own
nominal calibration and the falsifier did NOT fire:

    z = (D_native - expected)/sd      z_mean **-0.0521**   z_sd **1.9962**    (perfect: 0.00, 1.00)
    coverage of the central 50% band  **0.2824**   (nominal 0.50)
    coverage of the central 90% band  **0.6159**   (nominal 0.90)
    mean |D_native - expected| 2.3386 A against a published sd of 1.4471 A

**The posterior is very nearly CENTRED and roughly TWICE too narrow.**  Over-confidence peaks in the
mid-range separations (z_sd 1.23 at |i-j|=2, ~2.0 at 3-8, 1.28 at 9-15) and is flat in chain length.
24.1% of pairs carry a genuinely multimodal posterior, so any widening must preserve shape rather
than Gaussianise it.

WHY THIS SHOULD MOVE THE ENDPOINT, AND WHY IT IS NOT A REPEAT.  The shipped score is a Bayes risk
computed against the posterior, so an over-confident posterior yields an over-sharp risk curve and a
score that over-commits to a wrong distance.  Two independent prior results say that is expensive:
`confidently wrong costs 2-3x absent` from project memory, and s24 L13-A's measured DIRAC crossover
-- collapsing the posterior HURTS while the prior is wrong (+0.092 at gamma=0.1) and only helps once
it is nearly right (-0.049 at gamma=0.7), with the sign flipping at gamma ~ 0.55.  The shipped prior
sits at gamma = 0, deep in the region where sharpening is harmful.  **This arm moves in the opposite
direction from s24/conf.py, which reweighted PAIRS by confidence and failed; this reshapes each
pair's POSTERIOR.  Different object, opposite direction, mechanism stated in advance.**

THE ARMS.  Every arm is ONE GLOBAL PARAMETER -- the only router class the standing finite-sample
bound permits -- and every arm contains the incumbent as an exact interior point.

    SD(f)      per-pair distance-space convolution sized so the posterior's sd is multiplied by f.
               K[c,c'] proportional to exp(-(C[c']-C[c])^2 / 2 s^2) with s = sd_pair * sqrt(f^2-1),
               row-normalised, applied as p' = p @ K.  f = 1 is the identity.  PRIMARY: it targets the
               measured defect directly and adapts to each pair's own width.
    TEMP(T)    p' proportional to p^(1/T), renormalised.  T = 1 is the identity.  Widens by flattening
               relative weights; preserves mode LOCATIONS but not their spacing.  DECLARED SECONDARY.
    SDFIXW(f)  identical posterior to SD(f) but with the per-pair weight w frozen at the real prior's
               value, isolating "the posterior reshaped" from "the pair weighting moved".

FITTING, AND THE TEST THAT MAKES THIS MORE THAN A SWEEP.  Two nested-CV arms over the 5 pinned folds:

    CAL   the parameter is chosen on the TRAINING folds so that their z_sd is closest to 1.0.
          It never sees an RMSD.  If calibration is the mechanism, the calibration-optimal parameter
          must also improve held-out RMSD -- that is a genuine prediction, not a fit.
    RMSD  the parameter is chosen on the training folds by their mean RMSD, the usual construction.

**CAL is the primary.** A sweep that finds its best point by RMSD proves only that a grid has a
minimum; a parameter chosen by an independent criterion that then moves the endpoint is evidence
about the mechanism.

NATIVE USE.  The CAL arm uses TRAINING-FOLD native distances to estimate z_sd, and is applied to
held-out folds.  This is cross-target empirical calibration, the pattern the endgame directive SS10
names, and it is strictly less exposed than the predictor it corrects -- each of the 5 distogram fold
models already saw ~100 of the other 125 dev natives (s24 Workstream A).  No target's own native
informs its own correction.  Every held-out number is scored on targets whose natives were not used.

OPERATOR FORKS, per rule 0.  NOTE: enumerated by the coordinator, who has a stake; the audit lane
should re-enumerate.

    functional     DECLARED the shipped Bayes-risk functional, unmodified, recomputed from the widened
                   posterior by CONSTRUCTING A GENUINE `core.predict.Distogram` and using its own
                   `_risk` -- never a re-implementation.  Reconstruction from the shipped `prob` is
                   asserted bit-exact per target before any number is read.
                   NOT TAKEN re-implementing the risk table, and NOT TAKEN modifying the loss.
    basis          DECLARED point cloud throughout, medoid frame, as deployed.  NOT TAKEN a built chain.
    readout        DECLARED the shipped uniform top-75 coordinate average, so the ONLY thing that
                   varies is which 75 candidates are selected.  NOT TAKEN re-optimising m or weighting
                   members (s23 L5: uniform is optimal).
    normalisation  DECLARED the widening parameterised as a MULTIPLIER on each pair's own sd, so a
                   confident pair and a diffuse pair are widened proportionally.  NOT TAKEN a single
                   absolute width in Angstroms, which would erase the per-pair confidence structure
                   the posterior does carry.
    null           DECLARED f = 1 / T = 1, the shipped posterior, on identical data -- does widening
                   beat NOT widening.  NOT TAKEN the per-target oracle f, which is a best-of-K arm and
                   is reported separately, labelled ORACLE, and scored against its own order statistic.
    THE LABEL      DECLARED continuous Ca-RMSD, and SEPARATELY the held-out z_sd, reported together.
                   An RMSD gain WITHOUT z_sd moving toward 1 would refute the stated mechanism even if
                   the arm wins.  NOT TAKEN RMSD alone, and NOT TAKEN any binarised win rate.

  H  the posterior is over-confident, the score therefore over-commits, and widening it toward
     calibration improves held-out RMSD.
  Falsifier  the CAL arm failing to beat f=1 past its own MDE with a fold-clustered CI excluding zero.
  Secondary falsifier  an RMSD gain unaccompanied by z_sd moving toward 1.0, which would mean the arm
     wins for some reason other than the mechanism claimed.
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
from core import predict as dgm           # noqa: E402

TOPM = 75
FGRID = [1.0, 1.15, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0]
TGRID = [1.0, 1.25, 1.5, 2.0, 3.0, 5.0]


def _save(o, name="temper.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _widen_sd(prob, C, f):
    """Distance-space convolution sized per pair so the posterior's sd is multiplied by f."""
    if f <= 1.0:
        return prob.copy()
    m = (prob * C[None, :]).sum(1)
    v = (prob * (C[None, :] - m[:, None]) ** 2).sum(1)
    s = np.sqrt(np.maximum(v, 1e-12)) * np.sqrt(f * f - 1.0)
    d2 = (C[None, :] - C[:, None]) ** 2                       # (17,17)
    out = np.empty_like(prob)
    for p in range(len(prob)):
        K = np.exp(-d2 / (2.0 * max(s[p], 1e-6) ** 2))
        K /= K.sum(1, keepdims=True)
        q = prob[p] @ K
        out[p] = q / q.sum()
    return out


def _temper(prob, T):
    if T == 1.0:
        return prob.copy()
    q = np.power(np.maximum(prob, 1e-12), 1.0 / T)
    return q / q.sum(1, keepdims=True)


def _zsd(prob, C, Dt):
    m = (prob * C[None, :]).sum(1)
    v = (prob * (C[None, :] - m[:, None]) ** 2).sum(1)
    sd = np.sqrt(np.maximum(v, 1e-12))
    z = (Dt - m) / sd
    return float(z.std(ddof=1)) if len(z) > 1 else 0.0


def run():
    tg = I.targets()
    rows = []
    print("targets: %d  f grid %s  T grid %s" % (len(tg), FGRID, TGRID), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        seq = u["seq"]; nat = np.asarray(u["nat_ca"], float); n_res = int(u["n"])
        pi, pj = I.pair_index(n_res)
        dg = I.distogram(pdb, seq, u["fold"])
        prob0 = np.asarray(dg["prob"], float); C = np.asarray(dg["centres"], float)
        Wpool = np.asarray(u["W"], float)[I.pool_idx(u)]
        D = I.pair_dists(Wpool, pi, pj)
        Dt = I.pair_dists(nat[None], pi, pj)[0]

        #: GATE -- the shipped risk table must be reconstructible from the shipped posterior
        d0 = dgm.Distogram(seq, prob0, pi, pj)
        err = float(np.abs(np.asarray(dg["risk"], float) - np.asarray(d0._risk, float)).max())
        assert err < 1e-4, (pdb, err)
        w_real = np.asarray(d0.w, float)

        def emit(pr, fixw=False):
            d = dgm.Distogram(seq, pr, pi, pj)
            rk = np.asarray(d._risk, float)
            if fixw:
                rk = rk / np.maximum(np.asarray(d.w, float)[:, None], 1e-12) * w_real[:, None]
            sc = np.asarray(I.shipped_score({"grid": d.grid, "risk": rk}, D), float)
            sel = Wpool[np.argsort(sc, kind="stable")[:TOPM]]
            return float(I.ca_rmsd(_avg(sel), nat))

        row = {"pdb": pdb, "n": n_res, "fold": int(u["fold"]), "recon_err": err,
               "SD": [], "SDFIXW": [], "TEMP": [], "zsd_SD": [], "zsd_TEMP": []}
        for f in FGRID:
            pr = _widen_sd(prob0, C, f)
            row["SD"].append(emit(pr))
            row["SDFIXW"].append(emit(pr, fixw=True))
            row["zsd_SD"].append(_zsd(pr, C, Dt))
        for T in TGRID:
            pr = _temper(prob0, T)
            row["TEMP"].append(emit(pr))
            row["zsd_TEMP"].append(_zsd(pr, C, Dt))
        rows.append(row)
        if (c_i + 1) % 10 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("SD", "SDFIXW", "TEMP", "zsd_SD", "recon_err")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "fgrid": FGRID, "tgrid": TGRID, "topm": TOPM})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "temper.json")))["rows"]
    from s15 import seed as SD
    fold = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("temper", "rep")
    R = np.array([r["SD"] for r in rows], float)
    RF = np.array([r["SDFIXW"] for r in rows], float)
    RT = np.array([r["TEMP"] for r in rows], float)
    Z = np.array([r["zsd_SD"] for r in rows], float)
    ZT = np.array([r["zsd_TEMP"] for r in rows], float)
    base = R[:, 0]

    def st(x):
        x = np.asarray(x, float); n = len(x)
        se = x.std(ddof=1) / np.sqrt(n)
        F = sorted(set(fold))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se, float(np.percentile(fs, 2.5)),
                float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    print("\nn = %d.  Point cloud, shipped uniform top-75 average; only the SELECTION varies." % len(rows))
    print("  risk-table reconstruction from the shipped posterior: max err %.2e\n"
          % max(r["recon_err"] for r in rows))
    print("  %8s%11s%12s%11s%11s" % ("f", "RMSD", "RMSD fixw", "z_sd", "vs f=1"))
    for a, f in enumerate(FGRID):
        print("  %8.2f%11.4f%12.4f%11.4f%+11.4f"
              % (f, R[:, a].mean(), RF[:, a].mean(), Z[:, a].mean(), (R[:, a] - base).mean()))
    print("\n  %8s%11s%11s%11s" % ("T", "RMSD", "z_sd", "vs T=1"))
    for a, T in enumerate(TGRID):
        print("  %8.2f%11.4f%11.4f%+11.4f" % (T, RT[:, a].mean(), ZT[:, a].mean(),
                                              (RT[:, a] - RT[:, 0]).mean()))

    #: --- NESTED CV, two independent fitting criteria
    held_c = np.empty(len(rows)); held_r = np.empty(len(rows))
    fc = np.empty(len(rows)); fr = np.empty(len(rows))
    for q in sorted(set(fold)):
        tr = fold != q; te = fold == q
        a_cal = int(np.argmin(np.abs(Z[tr].mean(0) - 1.0)))     # never sees an RMSD
        a_rms = int(np.argmin(R[tr].mean(0)))
        held_c[te] = R[te][:, a_cal]; fc[te] = FGRID[a_cal]
        held_r[te] = R[te][:, a_rms]; fr[te] = FGRID[a_rms]

    print("\n  PRIMARY (paired vs the shipped posterior, f = 1):")
    for lab, x in (("CAL   f fitted on TRAINING-FOLD CALIBRATION  <-- PRIMARY", held_c - base),
                   ("RMSD  f fitted on training-fold RMSD", held_r - base),
                   ("per-target ORACLE f (best-of-%d, see null below)" % len(FGRID), R.min(1) - base)):
        m, se, mde, lo, hi, w = st(x)
        v = ("BEATS" if (hi < 0 and abs(m) > mde) else
             "worse" if (lo > 0 and abs(m) > mde) else "not measured")
        tm = " [TYPE-M]" if mde > 0 and 0.7 <= abs(m) / mde <= 1.3 else ""
        print("    %-52s%+.4f SE %.4f MDE %.4f fold[%+.4f,%+.4f] %3dW/%3dL %s%s"
              % (lab, m, se, mde, lo, hi, w, len(rows) - w, v, tm))
    print("    f selected by CALIBRATION: %s     by RMSD: %s"
          % (sorted(set(fc.tolist())), sorted(set(fr.tolist()))))
    print("    held-out mean RMSD, CAL arm: %.4f   (incumbent %.4f)" % (held_c.mean(), base.mean()))

    #: best-of-K null for the oracle row, within-target resampling
    dev = R - R.mean(1, keepdims=True)
    sim = np.array([dev[np.arange(len(rows)), rng.integers(0, len(FGRID), len(rows))]
                    for _ in range(2000)])
    nullgain = np.mean([np.minimum.reduce([dev[np.arange(len(rows)),
                        rng.integers(0, len(FGRID), len(rows))] for _ in range(len(FGRID))]).mean()
                        for _ in range(400)])
    obs = (R.min(1) - R.mean(1)).mean()
    print("    best-of-%d null on the ORACLE row: observed %+.4f, null %+.4f -> %.0f%% accounted"
          % (len(FGRID), obs, nullgain, 100 * nullgain / obs if obs else 0.0))

    print("\n  MECHANISM CHECK (the secondary falsifier):")
    zc = np.empty(len(rows))
    for q in sorted(set(fold)):
        tr = fold != q; te = fold == q
        zc[te] = Z[te][:, int(np.argmin(np.abs(Z[tr].mean(0) - 1.0)))]
    print("    held-out z_sd, shipped posterior %.4f  ->  CAL arm %.4f   (target 1.0)"
          % (Z[:, 0].mean(), zc.mean()))
    print("\n  Falsifier: the CAL arm failing to beat f=1 past its own MDE with a fold CI excluding zero.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
