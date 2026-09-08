"""TORSION-PREDICTOR -- extending Sprint 12's (sigma, coverage) surface. ORACLE DIAGNOSTIC.

Two things Sprint 12's surface cannot answer, both of which decide this agent's verdict:

1. **It stops at sigma = 20 degrees.** The published single-sequence state of the art sits at
   sigma ~ 41-44 degrees (lit §1.3, converted from MAE), and a leave-fold-out model on 79 000
   peptide/fragment residues will be worse still. Any read-off past 20 degrees is an
   extrapolation into a region where the emitted-RMSD axis is strongly convex, so the axis is
   measured here instead: sigma up to 100 degrees, same noise model, same builder, 3 seeds.

2. **Its missingness model is uniform random dropout.** The literature is specific that a real
   torsion channel's declined residues are CLUSTERED -- termini (a +/-3 heptapeptide matching
   window is truncated for 6 of 13 residues on a 13-mer) and contiguous dynamic stretches
   (RCI-S2 <= 0.6). A sequentially built chain does not care about the number of gaps, it
   cares about their run lengths: an isolated gap is bracketed by restrained neighbours, a
   contiguous run misplaces everything downstream. Three missingness models are measured at
   matched coverage -- `uniform`, `clustered` (geometric run lengths), `terminal` (gaps eaten
   from both ends first, then interior blocks).

EVERYTHING HERE READS NATIVE TORSIONS.  ORACLE DIAGNOSTIC, not a method, not deployable.

    python -m s13.tors_surface
"""
from __future__ import annotations
import os, sys, json, time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from s13 import tors_common as T                 # noqa: E402
from s13 import tors_eval as EV                  # noqa: E402
from s12 import instrument as I                  # noqa: E402

D2R = T.D2R
SIGMAS_FULL = (0, 6, 12, 20, 30, 40, 50, 60, 70, 80, 90, 100)
CELL_SIG = (12.0, 20.0, 30.0, 45.0)
CELL_FRAC = (1.0, 0.9, 0.75, 0.5)
SEEDS = (0, 1, 2)


# ------------------------------------------------------------------- missingness models
def mask_uniform(rng, n, frac):
    m = rng.random(n) < frac
    if not m.any():
        m[rng.integers(0, n)] = True
    return m


def mask_clustered(rng, n, frac, mean_run=3.0):
    """Drop `n*(1-frac)` residues as contiguous runs with geometric lengths."""
    ndrop = int(round((1 - frac) * n))
    m = np.ones(n, bool)
    guard = 0
    while (~m).sum() < ndrop and guard < 200:
        guard += 1
        run = max(1, int(rng.geometric(1.0 / mean_run)))
        run = min(run, ndrop - int((~m).sum()))
        s = int(rng.integers(0, max(1, n - run + 1)))
        m[s:s + run] = False
    if not m.any():
        m[rng.integers(0, n)] = True
    return m


def mask_terminal(rng, n, frac, term_share=0.7):
    """Gaps eaten from both termini first (the literature's dominant pattern), remainder as
    one interior block."""
    ndrop = int(round((1 - frac) * n))
    m = np.ones(n, bool)
    nt = int(round(term_share * ndrop))
    ln = int(rng.integers(0, nt + 1)) if nt > 0 else 0
    rn = nt - ln
    if ln:
        m[:ln] = False
    if rn:
        m[n - rn:] = False
    rem = ndrop - nt
    if rem > 0:
        lo = ln + 1
        hi = max(lo + 1, n - rn - rem)
        s = int(rng.integers(lo, hi)) if hi > lo else lo
        m[s:s + rem] = False
    if not m.any():
        m[rng.integers(0, n)] = True
    return m


MASKS = {"uniform": mask_uniform, "clustered": mask_clustered, "terminal": mask_terminal}


# ---------------------------------------------------------------------------- the sweeps
def run(nat):
    tg = I.targets()
    isf = np.array([t["pdb"] in set(I.FAIL18) for t in tg])
    out = {"full_coverage": {}, "cells": {}}

    for sg in SIGMAS_FULL:
        vals = []
        for t in tg:
            nphi, npsi = nat[t["pdb"]]
            _, _, ca = EV.pool_torsions(t["pdb"])
            acc = []
            for s in SEEDS:
                rng = np.random.default_rng(int(1e6 * sg + s + hash(t["pdb"]) % 9973))
                ph = T.wrap(nphi + rng.normal(0, sg * D2R, t["n"]))
                ps = T.wrap(npsi + rng.normal(0, sg * D2R, t["n"]))
                acc.append(I.ca_rmsd(I.build_ca(ph, ps), ca))
            vals.append(float(np.mean(acc)))
        v = np.array(vals)
        out["full_coverage"][f"{sg:g}"] = {
            "mean": float(v.mean()), "median": float(np.median(v)), "sd": float(v.std()),
            "frac_under_2.0": float((v < 2.0).mean()), "frac_under_1.5": float((v < 1.5).mean()),
            "FAIL18": float(v[isf].mean()), "other108": float(v[~isf].mean())}
        print(f"  sigma={sg:3g} full cov: {v.mean():.3f}  <2A {(v<2.0).mean():.2f}  "
              f"F18 {v[isf].mean():.3f}", flush=True)

    for sg in CELL_SIG:
        for fr in CELL_FRAC:
            for mk, fn in MASKS.items():
                if fr == 1.0 and mk != "uniform":
                    continue
                vals, runs = [], []
                for t in tg:
                    nphi, npsi = nat[t["pdb"]]
                    acc = []
                    for s in SEEDS:
                        rng = np.random.default_rng(int(1e6 * sg + 1e4 * fr * 100 + s
                                                        + hash(t["pdb"] + mk) % 9973))
                        ph = T.wrap(nphi + rng.normal(0, sg * D2R, t["n"]))
                        ps = T.wrap(npsi + rng.normal(0, sg * D2R, t["n"]))
                        m = fn(rng, t["n"], fr)
                        gaps = np.diff(np.flatnonzero(np.diff(np.r_[1, (~m).astype(int), 1])))
                        runs += [int(x) for x in gaps[::2] if x > 0]
                        pf, pp, ca = EV.fill_pool(t["pdb"], ph, ps, m)
                        acc.append(I.ca_rmsd(I.build_ca(pf, pp), ca))
                    vals.append(float(np.mean(acc)))
                v = np.array(vals)
                out["cells"][f"s{sg:g}_f{fr:g}_{mk}"] = {
                    "mean": float(v.mean()), "median": float(np.median(v)),
                    "frac_under_2.0": float((v < 2.0).mean()),
                    "FAIL18": float(v[isf].mean()), "other108": float(v[~isf].mean()),
                    "mean_gap_run": float(np.mean(runs)) if runs else 0.0}
                print(f"  sigma={sg:g} cov={fr:g} {mk:9s}: {v.mean():.3f} "
                      f"<2A {(v<2.0).mean():.2f} run={np.mean(runs) if runs else 0:.1f}",
                      flush=True)
    return out


def main():
    T.free_ok(1.5, tag="tors_surface")
    nat = EV.native_torsions()
    t0 = time.time()
    out = run(nat)
    out["what"] = ("ORACLE DIAGNOSTIC: Sprint 12's (sigma, coverage) surface extended to "
                   "sigma=100 deg and to CLUSTERED / TERMINAL missingness")
    out["secs"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(T.RESULTS, "tors_surface.json"), "w"), indent=1, default=str)
    print("wrote s13/results/tors_surface.json", f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
