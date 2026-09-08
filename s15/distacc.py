"""SPRINT 15, coordinator -- how accurate must the distances be, and is the sd calibrated?

CONTEXT. `s15/distgeo.py` established, by a route this project has never taken, that fitting
continuous torsions to a distance matrix is a GENERATION method that sidesteps the selection
bottleneck entirely. On a first 8-target read:

    ORACLE, fit to the TRUE distance matrix      0.697 mean, 0.077 median, 88% under 2 A
    predicted distances, inverse-variance        2.792
    predicted distances, unweighted              3.102
    the retrieval start it was given             3.247

Two things follow immediately. **The representation is not the barrier** -- torsion space
reproduces a true distance matrix essentially exactly, which is a stronger statement than the
k=4 descent ceiling because it is continuous and solved rather than searched. And **the
per-pair uncertainty is worth 0.31 A**, the first time this project has used the distogram's
`sd` for anything.

So the entire gap 0.697 -> 2.792 is distogram error. This module prices that gap.

THREE QUESTIONS, in the order they decide things.

**A. What IS the distogram's error?** Magnitude, dependence on sequence separation, and
whether it is biased or merely noisy. The locality theorem says a separation-8 pair is a
14-qubit interaction, so long-separation pairs are also the expensive ones -- if the error
concentrates there, that is both the problem and the opportunity.

**B. Is the `sd` calibrated?** If the standardised residual `(dhat - d_true)/sd` is not
approximately N(0,1), the inverse-variance weights are wrong and a recalibration is free
accuracy. Nobody has ever checked. An over-confident sd on long-range pairs would actively
mis-weight the fit toward exactly the pairs it gets wrong.

**C. How accurate would the distances have to be?** The phase diagram: corrupt the TRUE
matrix at a controlled error level, fit, and measure. Contours at 3.0 / 2.5 / 2.0 A, with the
distogram's ACTUAL measured error marked on it. This is the distance-channel analogue of the
torsion sigma-coverage diagram and it is a headline paper figure: it converts "the distogram
is not good enough" into "the distogram needs to be X% better, and here is where".

Error models are realistic, not just i.i.d.: relative and absolute Gaussian, a
separation-dependent model matched to the measured profile, a biased model, and sparse gross
outliers. A method that only reaches 2 A under unrealistically clean restraints is not a 2 A
method, and the brief requires that be said out loud.

ORACLE DIAGNOSTIC for the corruption arms and the error measurement (both read the native to
price a channel). The distogram itself is leave-fold-out and native-free; measuring its error
post hoc is evaluation, not leakage.

Run:
    python -m s15.distacc
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD             # noqa: E402
from s15 import distgeo as D                 # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

SEPS = [(2, 3), (4, 5), (6, 7), (8, 10), (11, 15)]
SIGMAS = [0.0, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00]


# ---------------------------------------------------------------- A/B: error and calibration
def error_profile(targets=None):
    """Distogram error and sd calibration, resolved by sequence separation."""
    tg = targets if targets is not None else I.targets()
    rows, allz, alle, alls, allsep = [], [], [], [], []
    for t in tg:
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        i, j = I.pair_index(n)
        dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
        dg = I.distogram(pdb, seq, fold)
        dhat = np.asarray(dg["expected"], float)
        sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)
        err = dhat - dtrue
        z = err / sd
        sep = (j - i)
        rows.append({"pdb": pdb, "n": n, "fold": fold,
                     "mae": float(np.abs(err).mean()), "bias": float(err.mean()),
                     "rmse": float(np.sqrt((err ** 2).mean())),
                     "z_mean": float(z.mean()), "z_sd": float(z.std()),
                     "sd_mean": float(sd.mean())})
        allz.append(z); alle.append(err); alls.append(sd); allsep.append(sep)
    z = np.concatenate(allz); e = np.concatenate(alle)
    s = np.concatenate(alls); sp = np.concatenate(allsep)

    by_sep = []
    for lo, hi in SEPS:
        m = (sp >= lo) & (sp <= hi)
        if not m.any():
            continue
        by_sep.append({"sep": f"{lo}-{hi}", "n_pairs": int(m.sum()),
                       "mae": float(np.abs(e[m]).mean()),
                       "bias": float(e[m].mean()),
                       "rmse": float(np.sqrt((e[m] ** 2).mean())),
                       "sd_mean": float(s[m].mean()),
                       "z_sd": float(z[m].std()),
                       "calibration": float(z[m].std())})
    return {"per_target": rows, "by_separation": by_sep,
            "global": {"mae": float(np.abs(e).mean()), "bias": float(e.mean()),
                       "rmse": float(np.sqrt((e ** 2).mean())),
                       "z_mean": float(z.mean()), "z_sd": float(z.std()),
                       "n_pairs": int(len(e))}}


# ------------------------------------------------------------------ C: the phase diagram
def corrupt(dtrue, sep, mode, sigma, rng):
    """ORACLE. Corrupt a true distance matrix under a named realistic error model."""
    d = dtrue.copy()
    if mode == "abs_gauss":
        d = d + rng.normal(0, sigma, size=d.shape)
    elif mode == "rel_gauss":
        d = d * (1.0 + rng.normal(0, sigma, size=d.shape))
    elif mode == "sep_scaled":
        # error grows with sequence separation, matching the shape a distogram usually has
        scale = sigma * (sep / sep.mean())
        d = d + rng.normal(0, 1, size=d.shape) * scale
    elif mode == "biased":
        d = d + sigma + rng.normal(0, sigma * 0.5, size=d.shape)
    elif mode == "outliers":
        d = d + rng.normal(0, sigma * 0.3, size=d.shape)
        k = max(1, int(0.05 * len(d)))
        idx = rng.choice(len(d), size=k, replace=False)
        d[idx] += rng.normal(0, sigma * 6.0, size=k)
    else:
        raise ValueError(mode)
    return np.maximum(d, 2.0)


def phase_diagram(targets=None, modes=("abs_gauss", "sep_scaled", "outliers"),
                  sigmas=SIGMAS, n_start=4, n_rep=2):
    tg = targets if targets is not None else I.targets()
    cells = {}
    for mode in modes:
        for sg in sigmas:
            vals = []
            for t in tg:
                pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
                u = I.load_univ(pdb)
                nat = np.asarray(u["nat_ca"], float)
                i, j = I.pair_index(n)
                dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
                sep = (j - i).astype(float)
                rr = []
                for rep in range(n_rep if sg > 0 else 1):
                    rng = SD.stable_rng(pdb, mode, sg, rep)
                    d = corrupt(dtrue, sep, mode, sg, rng) if sg > 0 else dtrue
                    w = np.ones_like(d)
                    S = D.starts(pdb, seq, n, fold, n_start, rng)
                    best = None
                    for phi0, psi0, _tag in S:
                        phi, psi, f, _ = D.fit_distances(d, w, i, j, phi0, psi0)
                        if best is None or f < best[0]:
                            best = (f, I.ca_rmsd(I.build_ca(phi, psi), nat))
                    rr.append(best[1])
                vals.append(float(np.mean(rr)))
            v = np.asarray(vals)
            cells[f"{mode}_s{sg:g}"] = {"mode": mode, "sigma": sg,
                                        "mean": float(v.mean()),
                                        "median": float(np.median(v)),
                                        "frac_under_2": float((v < 2.0).mean()),
                                        "frac_under_2_5": float((v < 2.5).mean()),
                                        "per_target": vals}
            print(f"  {mode:<12} sigma {sg:>5.2f} A -> {v.mean():6.3f} A "
                  f"(<2.5 {float((v<2.5).mean()):.2f}, <2.0 {float((v<2.0).mean()):.2f})",
                  flush=True)
    return cells


def threshold(cells, mode, target_rmsd):
    xs = sorted({c["sigma"] for c in cells.values() if c["mode"] == mode})
    ys = [next(c["mean"] for c in cells.values()
               if c["mode"] == mode and c["sigma"] == x) for x in xs]
    if target_rmsd <= ys[0]:
        return 0.0
    if target_rmsd >= ys[-1]:
        return float(xs[-1])
    return float(np.interp(target_rmsd, ys, xs))


def run(n_targets=40):
    tg = I.targets()[:n_targets]
    print(f"PART A/B -- distogram error and sd calibration, {len(I.targets())} targets")
    ep = error_profile()
    g = ep["global"]
    print(f"  global: MAE {g['mae']:.3f} A  bias {g['bias']:+.3f}  RMSE {g['rmse']:.3f}  "
          f"z-mean {g['z_mean']:+.3f}  z-sd {g['z_sd']:.3f}  ({g['n_pairs']} pairs)")
    print(f"  {'separation':<12}{'n':>8}{'MAE':>8}{'bias':>8}{'RMSE':>8}"
          f"{'sd':>8}{'z-sd':>8}")
    for b in ep["by_separation"]:
        print(f"  {b['sep']:<12}{b['n_pairs']:>8}{b['mae']:>8.3f}{b['bias']:>+8.3f}"
              f"{b['rmse']:>8.3f}{b['sd_mean']:>8.3f}{b['z_sd']:>8.3f}")
    print("  (z-sd = 1.0 means the sd is calibrated; >1 over-confident, <1 conservative)")

    print(f"\nPART C -- phase diagram, {len(tg)} targets")
    cells = phase_diagram(targets=tg)

    thr = {}
    for mode in {c["mode"] for c in cells.values()}:
        thr[mode] = {"3.0": threshold(cells, mode, 3.0),
                     "2.5": threshold(cells, mode, 2.5),
                     "2.0": threshold(cells, mode, 2.0)}
    out = {"error_profile": ep, "cells": cells, "thresholds": thr,
           "n_targets_phase": len(tg)}
    with open(os.path.join(RESULTS, "distacc.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_distacc", out, n_expected=len(tg))

    print(f"\nDistance accuracy required (A), by error model:")
    print(f"  {'mode':<14}{'for 3.0 A':>12}{'for 2.5 A':>12}{'for 2.0 A':>12}")
    for m, d in thr.items():
        print(f"  {m:<14}{d['3.0']:>12.2f}{d['2.5']:>12.2f}{d['2.0']:>12.2f}")
    print(f"\nThe distogram's ACTUAL error is MAE {g['mae']:.3f} A / RMSE {g['rmse']:.3f} A.")
    return out


if __name__ == "__main__":
    run()
