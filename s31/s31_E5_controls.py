#!/usr/bin/env python
"""s31/s31_E5_controls.py -- S31 lane E: the magnitude and direction controls for the split.

The decomposition arms of `s31_E2_deltas.py` differ in MAGNITUDE as well as in direction
(`rms ||y_along|| = 3.12` against `||y_perp|| = 1.99` against `||y|| = 3.70`), so an asymmetry
between them could be magnitude alone.  Contract rule 7: a control must match the operator's own
space.  Two families, both **ORACLE / NOT DEPLOYABLE**:

  1. SHRINK CURVE.  `c * y` at four c.  A perfect correction at reduced amplitude.  This prices
     "how much of an arm's effect is just how big it is".
  2. ENERGY-MATCHED DIRECTION SPLIT.  A direction `u` built to capture the SAME fraction of y's
     energy as `mu` does in that target, but otherwise arbitrary:
         u = cos(t) * y/||y|| + sin(t) * g_perp,   cos^2(t) = ||y_along(mu)||^2 / ||y||^2
     Then split y against `u` exactly as against `mu`.  If removing the energy-matched along
     component buys what removing `y_along(mu)` buys, then `mu` is not special and only the
     energy fraction matters.  If it does not, `mu` is special.
  3. RANDOM-DIRECTION SPLIT.  The unmatched null, for scale: an isotropic direction.

    python s31/s31_E5_controls.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s31 import s31_E_lib as L                     # noqa: E402
from s31.s31_E2_deltas import p_along, p_perp      # noqa: E402

OUT_NPZ = os.path.join(HERE, "results", "s31_E5_deltas.npz")
OUT_JSON = os.path.join(HERE, "results", "s31_E5_deltas.json")
SEED = 310505


def main():
    d = L.load()
    y, mu, tgt = d["y"], d["mu"], d["tgt"]
    rng = np.random.default_rng(SEED)
    A = {"PROD": np.zeros_like(y)}
    for c in (0.25, 0.50, 0.75, 0.90):
        A["SHRINK_Y_%03d" % int(c * 100)] = c * y
    em_a = np.empty_like(y); em_p = np.empty_like(y)
    rd_a = np.empty_like(y); rd_p = np.empty_like(y)
    frac = []
    for _, m in L.per_target(y, tgt):
        yy, uu = y[m], mu[m]
        ny = np.linalg.norm(yy)
        f = (np.linalg.norm(p_along(yy, uu)) / ny) ** 2 if ny > 1e-12 else 0.0
        frac.append(float(f))
        g = rng.standard_normal(len(yy))
        gp = g - (g @ yy) / max(yy @ yy, 1e-12) * yy          # component of g orthogonal to y
        gp = gp / max(np.linalg.norm(gp), 1e-12)
        u = np.sqrt(f) * (yy / max(ny, 1e-12)) + np.sqrt(max(1.0 - f, 0.0)) * gp
        em_a[m] = p_along(yy, u); em_p[m] = p_perp(yy, u)
        r = rng.standard_normal(len(yy))
        rd_a[m] = p_along(yy, r); rd_p[m] = p_perp(yy, r)
    A["ENERGYMATCHED_ALONG"] = em_a; A["ENERGYMATCHED_PERP"] = em_p
    A["RANDDIR_ALONG"] = rd_a; A["RANDDIR_PERP"] = rd_p

    diag = {"seed": SEED, "status": "ORACLE / NOT DEPLOYABLE -- every arm reads d_nat",
            "energy_fraction_along_mu_per_target_mean": float(np.mean(frac)),
            "arms": {k: {"rms": float(np.sqrt((v ** 2).mean())),
                         "norm_ratio_to_y": float(np.linalg.norm(v) / np.linalg.norm(y))}
                     for k, v in A.items()}}
    tmp = OUT_NPZ[:-4] + ".%d.tmp.npz" % os.getpid()
    np.savez_compressed(tmp, off=d["off"], pdbs=d["pdbs"], names=np.array(sorted(A)),
                        **{"delta_" + k: v for k, v in A.items()})
    os.replace(tmp, OUT_NPZ)
    tj = OUT_JSON + ".%d.tmp" % os.getpid()
    with open(tj, "w") as fh:
        json.dump(diag, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tj, OUT_JSON)
    for k in sorted(A):
        print("%-22s rms %.4f  norm/||y|| %.4f" % (k, diag["arms"][k]["rms"],
                                                   diag["arms"][k]["norm_ratio_to_y"]))
    print("mean per-target energy fraction of y along mu: %.4f"
          % diag["energy_fraction_along_mu_per_target_mean"])
    print("wrote", OUT_NPZ)


if __name__ == "__main__":
    main()
