"""SPRINT 15, INFO, PART C.4 -- the mechanism: real torsion errors live in the RMSD-quiet
subspace, and that is why the i.i.d. phase surface is the wrong instrument for them.

`s15.info_errstruct` showed that destroying only the SIGNS of a real emitter's angular
errors -- keeping every |error| at every position exactly -- costs +0.88 A (incumbent) to
+1.60 A (ORACLE best pool window).  So the information is in the DIRECTION of the error
vector in torsion space, not in its size.  This module names the direction.

Linearise the CA trace about the native torsions:

    J  =  d(superposed CA coordinates) / d(theta),   theta = (phi_1..phi_n, psi_1..psi_n)

computed by central differences at 1 degree and projected out of the six rigid-body modes,
so `||J e|| / sqrt(n)` is the first-order CA-RMSD produced by an angular error `e`.  For a
RANDOM unit direction the expected damage is `s_rms = sqrt(trace(J'J) / 2n)`.  Define

    alignment(e)  =  ||J e|| / (||e|| * s_rms)

= the damage a channel's error direction does, relative to a random direction of the same
magnitude.  alignment = 1 is i.i.d.-equivalent; alignment < 1 means the channel's error is
concentrated in the quiet directions and the i.i.d. surface over-prices it by exactly that
factor.

    python -m s15.info_null
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s15 import info_lib as L            # noqa: E402
from s15.info_regime import channels, wrap   # noqa: E402
import peptide_db as pdb                  # noqa: E402

H = np.radians(1.0)


def jacobian(phi, psi):
    """(3n, 2n) derivative of the SUPERPOSED CA trace wrt the torsions, central differences.

    Superposition is what makes the derivative meaningful: the six rigid-body directions
    carry no RMSD, so each perturbed trace is Kabsch-superposed onto the unperturbed one
    before the difference is taken.
    """
    n = len(phi)
    base = I.build_ca(phi, psi)
    cols = []
    for k in range(2 * n):
        dp = np.zeros(n); ds = np.zeros(n)
        (dp if k < n else ds)[k % n] = H
        a = I.build_ca(phi + dp, psi + ds)
        b = I.build_ca(phi - dp, psi - ds)
        a = I.superpose_batch(a[None], base)[0]
        b = I.superpose_batch(b[None], base)[0]
        cols.append(((a - b) / (2 * H)).ravel())
    return np.column_stack(cols)


def main():
    rows = {}
    per_t = []
    for t in I.targets():
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        n = len(nphi)
        J = jacobian(nphi, npsi)
        s_rms = float(np.sqrt((J ** 2).sum() / (2 * n)))
        rec = {"pdb": p, "n": n, "s_rms": s_rms}
        for nm, (phi, psi) in channels(t).items():
            e = np.concatenate([wrap(np.asarray(phi, float) - nphi),
                                wrap(np.asarray(psi, float) - npsi)])
            ne = float(np.linalg.norm(e))
            if ne < 1e-12:
                continue
            al = float(np.linalg.norm(J @ e) / (ne * s_rms))
            rows.setdefault(nm, []).append(al)
            rec[nm] = al
        per_t.append(rec)

    out = {"per_target": per_t, "alignment": {}}
    print("PART C.4 -- ALIGNMENT OF REAL TORSION ERRORS WITH THE RMSD-QUIET SUBSPACE\n")
    print("alignment = damage per unit angular error, relative to a random direction.")
    print("1.000 = i.i.d.-equivalent;  < 1 = the i.i.d. surface OVER-prices the channel.\n")
    print(f"{'channel':<34}{'mean':>8}{'median':>9}{'sd':>8}{'CI95':>20}"
          f"{'frac<1':>9}")
    for nm, v in rows.items():
        v = np.asarray(v, float)
        rb = np.random.default_rng(3)
        b = np.asarray([v[rb.integers(0, len(v), len(v))].mean() for _ in range(4000)])
        out["alignment"][nm] = {"mean": float(v.mean()), "median": float(np.median(v)),
                                "sd": float(v.std()), "n": int(len(v)),
                                "ci95": [float(np.percentile(b, 2.5)),
                                         float(np.percentile(b, 97.5))],
                                "frac_below_1": float((v < 1).mean())}
        s = out["alignment"][nm]
        print(f"{nm:<34}{s['mean']:>8.3f}{s['median']:>9.3f}{s['sd']:>8.3f}"
              f"   [{s['ci95'][0]:.3f},{s['ci95'][1]:.3f}]{s['frac_below_1']:>9.3f}")

    # NULL: a random direction, and a sign-flipped copy of the real one, through the same J
    print("\nNULLS through the same Jacobians:")
    rng = np.random.default_rng(0)
    randal, flipal = [], []
    for t in I.targets():
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        n = len(nphi)
        J = jacobian(nphi, npsi)
        s_rms = float(np.sqrt((J ** 2).sum() / (2 * n)))
        for _ in range(8):
            e = rng.normal(0, 1, 2 * n)
            randal.append(np.linalg.norm(J @ e) / (np.linalg.norm(e) * s_rms))
        ch = channels(t)
        phi, psi = ch["ORACLE best pool window"]
        e = np.concatenate([wrap(np.asarray(phi, float) - nphi),
                            wrap(np.asarray(psi, float) - npsi)])
        for _ in range(8):
            f = e * rng.choice([-1.0, 1.0], 2 * n)
            flipal.append(np.linalg.norm(J @ f) / (np.linalg.norm(f) * s_rms))
    out["null_random_direction"] = {"mean": float(np.mean(randal)),
                                    "sd": float(np.std(randal))}
    out["null_signflipped_oracle_window"] = {"mean": float(np.mean(flipal)),
                                             "sd": float(np.std(flipal))}
    print(f"  random Gaussian direction          {np.mean(randal):.3f} "
          f"(sd {np.std(randal):.3f})   [must be ~1.000 by construction]")
    print(f"  ORACLE window error, signs flipped {np.mean(flipal):.3f} "
          f"(sd {np.std(flipal):.3f})")
    L.jwrite("info_null", out)
    return out


if __name__ == "__main__":
    main()
