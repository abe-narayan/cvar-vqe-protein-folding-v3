"""s18/objceil.py -- THE OBJECTIVE'S OWN CEILING.  Is the gap the distogram, or the functional form?

WHY THIS EXISTS.  Sprint 17 established that refining the coordinate average toward the deployed
distance objective **loses 0.561 A** while cutting the objective 71%: the objective's optimum is
worse than the structure the pipeline already builds.  Sprint 18 asks whether a degree-1
truncation fixes that.  But there is a prior question that bounds the whole branch and that
nobody has measured:

    Is the objective mis-aimed because the DISTOGRAM is wrong, or because the FUNCTIONAL FORM
    -- a weighted sum of squared distance residuals over Ca pairs -- is wrong?

These have opposite prescriptions.  If the form is fine and the distogram is the problem, the
programme should improve the predictor and every objective-engineering experiment is a
distraction.  If the form is also broken, then no distogram however good would rescue it, and
degree-1, degree-2 and Legacy-augmented variants are all rearranging the wrong object.

THE MEASUREMENT.  Refine the same starting structure toward the same functional form, varying
only the distances it is asked to satisfy:

    d_alpha = (1 - alpha) * dhat + alpha * d_true,   alpha in {0, 0.25, 0.5, 0.75, 1}

alpha = 0 is the deployed objective (Sprint 17 measured its endpoint at 3.610 A).  alpha = 1 is
the SAME objective with perfect distances -- the functional form's own ceiling.  The curve
between them is the sprint's requirement function: **how accurate must the distogram be for this
objective to reach 2.5 A?**

EVERY ARM WITH alpha > 0 IS AN ORACLE DIAGNOSTIC.  It reads native distances to construct the
target of an optimisation.  Nothing here is deployable and nothing here may be quoted as a
result.  It is a ceiling and a requirement, in the same class as Sprint 17's design equation.

CONTROLS, both mandatory.

    shuffled   dhat's residual re-randomised in DIRECTION while its per-pair MAGNITUDE
               distribution is preserved.  This separates "the distogram's errors are large"
               from "the distogram's errors point somewhere specific".  Sprint 15 recorded the
               prior as "full amplitude, wrong direction"; if that is right, the shuffled arm
               should be no worse than alpha = 0, because a wrong direction is a wrong
               direction.  If shuffling is much WORSE, the deployed distogram's error direction
               is carrying real information after all.
    isotropic  dhat perturbed by isotropic noise matched to the true residual's RMS, which is
               the "same magnitude, no structure" reference.

PRE-REGISTRATION, written before the run.

  HYPOTHESIS.  alpha = 1 reaches well below 2.5 A -- the functional form is sound and the
  distogram is the whole problem.

  EXPECTED.  The programme's record says perfect distance knowledge caps near 1.95-2.0 A through
  this geometry, and that oracle distances give 0.36 A in-pool.  So I expect alpha = 1 to land
  somewhere in 0.5-2.0 A.  If it lands near 3.6 A the functional form itself is broken and that
  is the most important negative result available this sprint.

  SUCCESS.  A monotone curve with alpha = 1 materially below 2.5 A, target as the unit, paired
  fold-clustered interval excluding zero against alpha = 0.

  FALSIFIER.  If alpha = 1 does NOT reach materially below the coordinate average (3.048 A),
  **the entire objective branch is closed** -- degree-1 included -- because even a perfect
  distogram cannot make this functional form point at the native, and the sprint should say so.
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
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import seed as SD                   # noqa: E402

ALPHAS = (0.0, 0.25, 0.5, 0.75, 1.0)


def _rmsd_of(phi, psi, nat):
    return float(I.ca_rmsd(I.build_ca(phi, psi), nat))


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows = []
    t0 = time.time()
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
        #: ORACLE -- the native's own Ca-Ca distances on the same pair index
        dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
        rng = SD.stable_rng(pdb, "s18objceil")

        W = np.asarray(AV.top75_windows(pdb)[0], float)
        P = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

        e = {"avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
             "proj": _rmsd_of(phi0, psi0, nat),
             "resid_rms": float(np.sqrt(((dhat - dtrue) ** 2).mean())),
             "resid_mae": float(np.abs(dhat - dtrue).mean())}

        #: the alpha ladder -- same functional form, better distances
        for a in ALPHAS:
            da = (1.0 - a) * dhat + a * dtrue
            p_, q_, _f = A.fit(da, sd, i, j, phi0, psi0)
            e[f"a{a}"] = _rmsd_of(p_, q_, nat)

        #: CONTROL 1 -- preserve the residual MAGNITUDES, destroy their DIRECTION.
        #: If the record's "full amplitude, wrong direction" reading is right, this should be
        #: no worse than alpha = 0.
        r = dhat - dtrue
        sh = dtrue + rng.permutation(r) * rng.choice([-1.0, 1.0], size=len(r))
        sh = np.maximum(sh, 2.0)
        p_, q_, _f = A.fit(sh, sd, i, j, phi0, psi0)
        e["shuffled"] = _rmsd_of(p_, q_, nat)

        #: CONTROL 2 -- isotropic noise matched to the true residual's RMS
        iso = np.maximum(dtrue + rng.standard_normal(len(r)) * e["resid_rms"], 2.0)
        p_, q_, _f = A.fit(iso, sd, i, j, phi0, psi0)
        e["isotropic"] = _rmsd_of(p_, q_, nat)

        #: CONTROL 3 -- SEPARATION-STRATIFIED shuffle.  The plain shuffle permutes residuals
        #: across all pairs, which also destroys their correlation with sequence separation --
        #: and the programme's record says the distogram's bias GROWS with separation.  So a
        #: plain shuffle may win by accidentally shrinking the long-range residuals the
        #: objective weights heavily, rather than by destroying "direction".  Permuting only
        #: WITHIN separation bins removes that confound.  If the effect dies here, the
        #: coordinator's headline is wrong.
        sep = np.abs(np.asarray(j) - np.asarray(i))
        r2 = r.copy()
        for lo, hi in ((2, 3), (4, 5), (6, 7), (8, 10), (11, 99)):
            m = (sep >= lo) & (sep <= hi)
            if m.sum() > 1:
                r2[m] = rng.permutation(r[m]) * rng.choice([-1.0, 1.0], size=int(m.sum()))
        st = np.maximum(dtrue + r2, 2.0)
        p_, q_, _f = A.fit(st, sd, i, j, phi0, psi0)
        e["shuf_strat"] = _rmsd_of(p_, q_, nat)

        #: CONTROL 4 -- shuffle the residual AND its weight together, so any adverse alignment
        #: between large residuals and large 1/sd^2 weights is preserved rather than broken.
        pi = rng.permutation(len(r))
        pw = np.maximum(dtrue + r[pi] * rng.choice([-1.0, 1.0], size=len(r)), 2.0)
        p_, q_, _f = A.fit(pw, sd[pi], i, j, phi0, psi0)
        e["shuf_paired"] = _rmsd_of(p_, q_, nat)

        rows.append({"pdb": pdb, "n": n, "fold": fold, **e})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": False},
                      open(os.path.join(RESULTS, "objceil.json"), "w"))
    json.dump({"rows": rows, "complete": len(rows) == len(tg)},
              open(os.path.join(RESULTS, "objceil.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "objceil.json")))["rows"]
    rng = SD.stable_rng("objceil", "report")
    g = lambda k: np.array([r[k] for r in rows])      # noqa: E731
    avg = g("avg"); a0 = g("a0.0")
    print(f"\nn = {len(rows)}.  Start = the coordinate average of the shipped top-75.")
    print("EVERY alpha > 0 ARM IS AN ORACLE DIAGNOSTIC: it builds the optimisation target from")
    print("native distances.  These are CEILINGS and REQUIREMENTS, never results.\n")
    print(f"  coordinate average (start)   {avg.mean():.3f}   median {np.median(avg):.3f}")
    print(f"  ideal-geometry projection    {g('proj').mean():.3f}")
    print(f"  distogram residual RMS       {g('resid_rms').mean():.3f} A   "
          f"MAE {g('resid_mae').mean():.3f} A\n")
    print(f"  {'arm':<14}{'RMSD':>8}{'median':>9}{'vs alpha=0':>24}{'vs the average':>24}{'W/L vs avg':>12}")
    order = [f"a{a}" for a in ALPHAS] + ["shuffled", "shuf_strat", "shuf_paired",
                                          "isotropic"]
    for k in order:
        v = g(k)
        m1, l1, h1 = _boot(v - a0, rng)
        m2, l2, h2 = _boot(v - avg, rng)
        w = int((v < avg).sum()); l = int((v > avg).sum())
        lab = k if not k.startswith("a") else f"alpha = {k[1:]}"
        tag = "  CONTROL" if k.startswith(("shuf", "isot")) else ""
        print(f"  {lab:<14}{v.mean():>8.3f}{np.median(v):>9.3f}"
              f"   {m1:+.3f} [{l1:+.3f},{h1:+.3f}]   {m2:+.3f} [{l2:+.3f},{h2:+.3f}]"
              f"{w:>6}/{l}{tag}")

    a1 = g("a1.0")
    print(f"\n  THE FUNCTIONAL FORM'S OWN CEILING (alpha = 1, perfect distances): "
          f"{a1.mean():.3f} A")
    print(f"  fraction of targets below 2.5 A at alpha = 1: {(a1 < 2.5).mean():.1%}")
    print(f"  fraction below 2.0 A: {(a1 < 2.0).mean():.1%}")
    print("\nREAD.  If alpha = 1 lands well below 2.5 A the functional form is sound and the whole")
    print("gap is DISTOGRAM ERROR -- the programme should improve the predictor and every")
    print("objective-engineering variant is a distraction.  If alpha = 1 fails to beat the")
    print("coordinate average, the FUNCTIONAL FORM itself is broken, no distogram would rescue")
    print("it, and the entire objective branch -- degree-1 included -- is closed.")
    print("\nThe `shuffled` control decides whether the deployed distogram's error DIRECTION")
    print("carries information: same magnitudes, randomised signs and pairing.  If shuffling is")
    print("no worse than alpha = 0, the record's 'full amplitude, wrong direction' reading holds.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
