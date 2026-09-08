"""SPRINT 19, AGENT A -- the shared fit harness, byte-identical to `s18/objceil.py`'s.

Start  = ideal-geometry projection of the coordinate average of the shipped top-75.
Fit    = `s15/align_lib.fit(d, sd, i, j, phi0, psi0)`, i.e. the DEPLOYED objective
         `sum w_ij (d_ij(phi,psi) - d_ij^target)^2` with `w = 1/sd^2`.
Metric = full-chain Ca-RMSD against model 1 of the native.

The start is deterministic and expensive (2-4 s) and identical for every arm, so it is
cached once per target.  `verify()` reproduces `objceil`'s `avg`, `proj` and `a0.0` columns.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import align_lib as A               # noqa: E402

CACHE = os.path.join(HERE, "cache")
os.makedirs(CACHE, exist_ok=True)


def start(pdb, seq, fold):
    """(phi0, psi0, avg_rmsd_inputs) -- cached ideal-geometry projection of the top-75 average."""
    path = os.path.join(CACHE, f"start_{pdb}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return z["phi"], z["psi"], z["avg"]
    W = np.asarray(AV.top75_windows(pdb)[0], float)
    P = I.pairwise_rmsd(W)
    avg, _b = I.coordinate_average(W, P)
    pr = I.project(np.asarray(avg, float), seq, int(fold))
    phi = np.asarray(pr["phi"], float)
    psi = np.asarray(pr["psi"], float)
    np.savez(path, phi=phi, psi=psi, avg=np.asarray(avg, float))
    return phi, psi, np.asarray(avg, float)


def fit_rmsd(target_d, sd, i, j, phi0, psi0, nat, floor=2.0):
    """One arm: fit to `target_d` with weights 1/sd^2 and return (rmsd, final objective)."""
    t = np.maximum(np.asarray(target_d, float), floor)
    ph, ps, f = A.fit(t, sd, i, j, phi0, psi0)
    return float(I.ca_rmsd(I.build_ca(ph, ps), nat)), float(f)


def obj_at(target_d, sd, i, j, d_of_X):
    """The deployed objective evaluated at a given distance vector."""
    w = 1.0 / np.asarray(sd, float) ** 2
    return float((w * (np.asarray(d_of_X, float) - np.asarray(target_d, float)) ** 2).sum())
