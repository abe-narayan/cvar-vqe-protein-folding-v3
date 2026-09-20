#!/usr/bin/env python
"""s29/s29_P_bandsel.py -- DIAGNOSTIC: is separation-band re-weighting of the OBJECTIVE the same
operation as rescaling the emitted cloud?

The coordinator's question (2026-09-20 00:45) after S29-L22 result 5: lane T offered
"separation-band re-weighting" as the only live form of the calibration item, and asked whether
re-weighting the objective's pairs BEFORE selection is genuinely different from rescaling the
emitted cloud AFTER it, or the same operation seen twice.

THE TEST, and its falsifier, registered here before the numbers:
Production scores a candidate by the UNIFORM mean over pairs of the distogram's per-pair L1 Bayes
risk (`s12.instrument.shipped_score`). Replace that uniform mean with a weighted mean whose
weights depend only on the sequence separation |i - j|, re-select the top 75 by the re-weighted
score (same tie key, same m), coordinate-average, and measure the emitted cloud's SEPARATION
PROFILE (S29-L22 result 4) and its RMSD.

  FALSIFIED (the operators are DIFFERENT, and result 5 does not reach T's route) if any band
  weighting moves the emitted profile's short-range end (the |i-j| = 1 ratio) by more than 0.02
  from production's 0.773, i.e. if selection can materially reshape the distortion.
  CONFIRMED (the profile is a property of AVERAGING, not of SELECTION) if every weighting -- including
  ones that put all the weight on the local band or all of it on the long-range band -- leaves the
  profile's shape essentially where production puts it.

This is a DIAGNOSTIC: no deployable parameter is chosen here and the built chain is not run.
ORACLE: the emitted cloud's RMSD reads the native, post hoc, for reporting only.

    python s29/s29_P_bandsel.py [--limit N]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s29 import s29_P_scale as P           # noqa: E402

OUT = os.path.join(P.RESULTS, "s29_P_bandsel.json")

#: weight(sep) schemes. "uniform" IS production and must reproduce it exactly.
SCHEMES = {
    "uniform": lambda sep: np.ones_like(sep, float),
    "local_only": lambda sep: (sep <= 3).astype(float),
    "long_only": lambda sep: (sep >= 6).astype(float),
    "local_x4": lambda sep: np.where(sep <= 3, 4.0, 1.0),
    "long_x4": lambda sep: np.where(sep >= 6, 4.0, 1.0),
    "inv_sep": lambda sep: 1.0 / sep.astype(float),
    "sep": lambda sep: sep.astype(float),
}


def weighted_score(dg, D, w):
    """The shipped score with a weighted mean over pairs in place of the uniform one.
    `w` is (npairs,) and non-negative; w == 1 reproduces `I.shipped_score` exactly."""
    grid = np.asarray(dg["grid"], float)
    risk = np.asarray(dg["risk"], float)
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    r = risk[np.arange(risk.shape[0])[None, :], g]          # (b, npairs)
    w = np.asarray(w, float)
    return (r * w[None, :]).sum(1) / w.sum()


def main(limit=None):
    fac = json.load(open(P.FACTORS, encoding="utf-8"))
    order = fac["order"][:limit] if limit else fac["order"]
    acc = {k: {"num": {}, "den": {}, "rmsd": [], "bond": [], "rg": [],
               "overlap": []} for k in SCHEMES}
    prod_check = []
    for pdb in order:
        cand, ch, _ = RP.channels_for(pdb)
        dg = I.distogram(pdb, cand.seq, cand.fold)
        i, j = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
        sep = np.abs(j - i)
        W = np.asarray(cand.W, float)
        D = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=-1)        # (k, npairs)
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        top_prod = np.lexsort((key, RP.zr(ch["DIS"])))[:RP.M]
        n = int(cand.n)
        nat = np.asarray(cand.nat_ca, float)                         # ORACLE, scoring only
        Dt = np.linalg.norm(nat[:, None, :] - nat[None, :, :], axis=-1)
        for name, fn in SCHEMES.items():
            w = fn(sep)
            E = weighted_score(dg, D, w)
            top = np.lexsort((key, RP.zr(E)))[:RP.M]
            C, _ = I.coordinate_average(W[top])
            if name == "uniform":
                prod_check.append(abs(float(I.ca_rmsd(C, nat))
                                      - float(fac["per_target"][pdb]["rmsd_cloud"])))
            a = acc[name]
            a["rmsd"].append(float(I.ca_rmsd(C, nat)))
            a["bond"].append(P.mean_bond(C)); a["rg"].append(P.rg(C))
            a["overlap"].append(len(set(top.tolist()) & set(top_prod.tolist())) / RP.M)
            Dc = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=-1)
            for s in range(1, n):
                aa = np.arange(n - s); bb = aa + s
                a["num"][s] = a["num"].get(s, 0.0) + float(Dc[aa, bb].sum())
                a["den"][s] = a["den"].get(s, 0.0) + float(Dt[aa, bb].sum())
    print("uniform arm reproduces production's cloud RMSD: max |diff| %.3e on %d targets"
          % (max(prod_check), len(prod_check)))
    print()
    seps = sorted(acc["uniform"]["num"])
    seps = [s for s in seps if acc["uniform"]["den"][s] > 0][:13]
    print("EMITTED SEPARATION PROFILE (cloud/native, ratio of means) BY SELECTION SCHEME")
    print("  %-11s %7s %7s %7s | %s" % ("scheme", "RMSD", "bond", "top75", "  ".join("%5d" % s for s in seps)))
    rows = {}
    for name in SCHEMES:
        a = acc[name]
        prof = [a["num"][s] / a["den"][s] for s in seps]
        rows[name] = dict(rmsd=float(np.mean(a["rmsd"])), bond=float(np.mean(a["bond"])),
                          rg=float(np.mean(a["rg"])), overlap=float(np.mean(a["overlap"])),
                          profile=prof, seps=seps, per_target_rmsd=a["rmsd"])
        print("  %-11s %7.4f %7.4f %7.3f | %s"
              % (name, rows[name]["rmsd"], rows[name]["bond"], rows[name]["overlap"],
                 "  ".join("%5.3f" % v for v in prof)))
    base = rows["uniform"]["profile"][0]
    worst = max(abs(rows[k]["profile"][0] - base) for k in SCHEMES)
    print()
    print("production's |i-j|=1 ratio %.4f ; largest move under ANY band weighting %.4f"
          % (base, worst))
    print("FALSIFIER (>0.02 move = the operators differ): %s"
          % ("FIRES -- selection CAN reshape the distortion" if worst > 0.02 else
             "DOES NOT FIRE -- the profile is a property of AVERAGING, not of SELECTION"))
    folds = np.array([fac["per_target"][p]["fold"] for p in order])
    print()
    for name in SCHEMES:
        if name == "uniform":
            continue
        o = ST.compare(np.array(rows[name]["per_target_rmsd"]),
                       np.array(rows["uniform"]["per_target_rmsd"]), folds=folds, names=order,
                       label="P band-select %s - production (POINT CLOUD, diagnostic)" % name)
        print(ST.fmt(o)); print()
    ST.save_atomic(OUT, dict(kind="separation-band re-weighting of the selection objective",
                             order=order, rows=rows, base_ratio=base, worst_move=worst,
                             falsifier_fires=bool(worst > 0.02),
                             uniform_reproduces_max_diff=float(max(prod_check))),
                   module_file=__file__)
    print("wrote", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--limit", type=int, default=None)
    main(**vars(ap.parse_args()))
