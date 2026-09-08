"""SPRINT 15, coordinator -- the distogram is BIASED, and the bias is correctable for free.

WHAT `s15/distacc.py` MEASURED, on all 126 targets and 8,549 pairs:

    separation      n     MAE    bias    RMSE    z-sd
    2-3          2636   0.974  +0.113   1.463   2.186
    4-5          2132   2.148  +0.323   2.963   2.983
    6-7          1628   2.818  +0.615   4.001   2.699
    8-10         1506   3.745  +0.931   5.216   2.701
    11-15         647   4.666  +1.492   6.328   2.671
    GLOBAL       8549   2.386  +0.509   3.704   2.633

Two defects, and only one of them matters for a weighted fit.

**The `sd` is 2.633x over-confident** -- the standardised residual has sd 2.633 where a
calibrated predictor gives 1.0. That is a real reporting defect and it must be stated in any
paper that uses this distogram. But for the FIT it is nearly harmless: the z-sd is roughly
CONSTANT across separation (2.19, 2.98, 2.70, 2.70, 2.67), so the sd captures the error's
SHAPE correctly and is only wrong by a near-uniform factor -- and a weighted least-squares
argmin is invariant under uniform rescaling of all weights. Recalibration alone changes
nothing about the structure. That is worth knowing before anyone spends a sprint on it.

**The BIAS does matter, and nobody has ever corrected it.** The distogram systematically
OVER-predicts distance, and the over-prediction grows monotonically with sequence separation,
from +0.113 A at separation 2-3 to **+1.492 A at separation 11-15**. A least-squares fit to
systematically inflated distances produces a systematically over-extended structure. This is
a pure, removable, additive error.

It is also nicely opposed to a known Sprint 14 result: coordinate averaging *contracts* the
backbone by 25.8%. The two stages have opposite geometric pathologies, which is worth stating.

**Correcting it is legitimate and native-free at inference.** The bias is a property of the
PREDICTOR, estimable from training folds alone, so a leave-fold-out bias profile `b(sep)`
applied as `d_corrected = dhat - b(sep)` uses no information about the target being predicted.
The five folds are already pinned. This module fits `b` on the four training folds and
applies it to the held-out fold, exactly as every other learned component in this project.

ARMS
  raw                  the distogram as shipped, inverse-variance weighted (the s15/distgeo
                       predictive arm, for reference)
  debias_global        one scalar correction, leave-fold-out
  debias_sep           a separation-resolved correction, leave-fold-out  <- the real arm
  debias_sep_affine    a separation-resolved affine map d -> a(sep)*d + b(sep), leave-fold-out
  ORACLE_debias_sep    the same correction fitted ON the held-out fold -- the ceiling of
                       debiasing, and an ORACLE DIAGNOSTIC that can never be a headline

The ORACLE arm is the control that says whether the leave-fold-out correction is limited by
transfer or by the correction's functional form -- a distinction Sprint 14 showed matters
enormously (in-band ordering was learnable at 0.986 within a target and 0.600 across).

Run:
    python -m s15.distcal
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
CACHE = os.path.join(ROOT, "s15", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

SEP_EDGES = [2, 4, 6, 8, 11, 999]          # bins: 2-3, 4-5, 6-7, 8-10, 11+


def sep_bin(sep):
    return np.clip(np.searchsorted(SEP_EDGES, sep, side="right") - 1, 0, len(SEP_EDGES) - 2)


def gather(targets):
    """Per-target predicted / true distances and separations. ORACLE for `dtrue` only."""
    out = {}
    for t in targets:
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        i, j = I.pair_index(n)
        dg = I.distogram(pdb, seq, fold)
        out[pdb] = {
            "i": i, "j": j, "n": n, "fold": fold, "seq": seq,
            "dhat": np.asarray(dg["expected"], float),
            "sd": np.maximum(np.asarray(dg["sd"], float), 1e-3),
            "dtrue": np.sqrt(((nat[i] - nat[j]) ** 2).sum(1)),
            "sep": (j - i).astype(float), "nat": nat,
        }
    return out


def fit_correction(data, pdbs, kind):
    """Fit a bias correction on the given targets only. Returns a callable."""
    dh = np.concatenate([data[p]["dhat"] for p in pdbs])
    dt = np.concatenate([data[p]["dtrue"] for p in pdbs])
    sp = np.concatenate([data[p]["sep"] for p in pdbs])
    if kind == "global":
        b = float((dh - dt).mean())
        return lambda d, s: d - b, {"b": b}
    if kind == "sep":
        bins = sep_bin(sp)
        b = np.array([float((dh[bins == k] - dt[bins == k]).mean())
                      if (bins == k).any() else 0.0
                      for k in range(len(SEP_EDGES) - 1)])
        return (lambda d, s: d - b[sep_bin(s)]), {"b_by_bin": b.tolist()}
    if kind == "sep_affine":
        bins = sep_bin(sp)
        A = np.ones((len(SEP_EDGES) - 1, 2))
        for k in range(len(SEP_EDGES) - 1):
            m = bins == k
            if m.sum() > 10:
                # least squares  dtrue ~ a*dhat + c
                X = np.column_stack([dh[m], np.ones(m.sum())])
                coef, *_ = np.linalg.lstsq(X, dt[m], rcond=None)
                A[k] = coef
        return (lambda d, s: A[sep_bin(s), 0] * d + A[sep_bin(s), 1]), {"affine": A.tolist()}
    raise ValueError(kind)


def run(targets=None, n_start=6):
    tg = targets if targets is not None else I.targets()
    data = gather(tg)
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([data[p]["fold"] for p in pdbs], int)
    fail = np.isin(pdbs, I.FAIL18)

    # ---- leave-fold-out corrections, and the ORACLE in-fold ceiling
    corr = {}
    for kind in ("global", "sep", "sep_affine"):
        corr[kind] = {}
        for f in sorted(set(folds)):
            train = [p for p in pdbs if data[p]["fold"] != f]
            corr[kind][f], _ = fit_correction(data, train, kind)
    orc = {}
    for f in sorted(set(folds)):
        held = [p for p in pdbs if data[p]["fold"] == f]
        orc[f], _ = fit_correction(data, held, "sep")

    arms = ["raw", "debias_global", "debias_sep", "debias_sep_affine", "ORACLE_debias_sep"]
    res = {a: [] for a in arms}
    path = os.path.join(RESULTS, "distcal.json")

    for c, p in enumerate(pdbs):
        d = data[p]
        rng = SD.stable_rng(p)
        S = D.starts(p, d["seq"], d["n"], d["fold"], n_start, rng)
        w = 1.0 / (d["sd"] ** 2)
        f = d["fold"]
        targets_by_arm = {
            "raw": d["dhat"],
            "debias_global": corr["global"][f](d["dhat"], d["sep"]),
            "debias_sep": corr["sep"][f](d["dhat"], d["sep"]),
            "debias_sep_affine": corr["sep_affine"][f](d["dhat"], d["sep"]),
            "ORACLE_debias_sep": orc[f](d["dhat"], d["sep"]),
        }
        for a in arms:
            tgt = np.maximum(targets_by_arm[a], 2.0)
            best = None
            for phi0, psi0, _t in S:
                phi, psi, fv, _ = D.fit_distances(tgt, w, d["i"], d["j"], phi0, psi0)
                if best is None or fv < best[0]:
                    best = (fv, float(I.ca_rmsd(I.build_ca(phi, psi), d["nat"])))
            res[a].append(best[1])
        if (c + 1) % 10 == 0 or c + 1 == len(pdbs):
            with open(path, "w") as fh:
                json.dump({"partial": {a: res[a] for a in arms}, "n_done": c + 1,
                           "pdbs": pdbs[:c + 1]}, fh)
            print(f"  {c+1}/{len(pdbs)} checkpointed", flush=True)

    from s14 import ladder as L
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    raw = np.asarray(res["raw"], float)

    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "arms": {}}
    for a in arms:
        v = np.asarray(res[a], float)
        out["arms"][a] = {
            **I.summary(v), "median": float(np.median(v)),
            "FAIL18": float(v[fail].mean()),
            "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs),
            "vs_raw": I.paired(v, raw, folds=folds, names=pdbs),
        }
    out["per_target"] = {a: dict(zip(pdbs, map(float, res[a]))) for a in arms}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_distcal", out, n_expected=len(pdbs))

    print(f"\nincumbent {ref.mean():.3f}\n")
    print(f"{'arm':<22}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs raw':>24}{'vs incumbent':>24}")
    for a in arms:
        s = out["arms"][a]
        vr, vi = s["vs_raw"], s["vs_incumbent"]
        u25 = float((np.asarray(res[a]) < 2.5).mean())
        print(f"{a:<22}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{u25:>8.2f}{s['FAIL18']:>9.3f}"
              f"  {vr['mean_diff']:+.3f} [{vr['ci95'][0]:+.3f},{vr['ci95'][1]:+.3f}]"
              f"  {vi['mean_diff']:+.3f} [{vi['ci95'][0]:+.3f},{vi['ci95'][1]:+.3f}]")
    return out


if __name__ == "__main__":
    run()
