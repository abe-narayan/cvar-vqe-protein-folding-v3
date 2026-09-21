"""LANE L / L1b -- the representability floor is LENGTH, not corpus.

`s32_L_feasibility census` measured the ideal-geometry rebuild error at 0.347 A on the
126 nine-to-sixteen-residue peptides and 2.53-4.66 A on 46-80 residue proteins.  Those are
two different CORPORA as well as two different LENGTHS, so the census alone cannot say
which one produced the rise.  The project's most repeated error is a control that does not
match the operator's space (memory: control-must-match-the-operators-space).

The control here is exact.  For ONE protein, take contiguous sub-windows of length L for
several L, rebuild each from its OWN native torsions, and measure the CA-RMSD.  Molecule,
resolution, corpus, deposition year and admission policy are all held fixed; only L moves.
If the floor rises with L inside a single molecule, the rise is length.

SELF-TEST (contract rule 5 -- a verification must be able to fail).  `selftest()` feeds the
same estimator a chain BUILT BY the ideal-geometry builder from random torsions.  Such a
chain lies exactly on the manifold, so a correct estimator must return ~0 at EVERY length.
An estimator that simply grew with L -- the failure mode that would fake this result --
returns a rising curve here and fails the assert.

    python -m s32.s32_L_floor run
    python -m s32.s32_L_floor selftest
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

from core import geometry as geo                                            # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROT_DIR = os.path.join(BASE, "prots")
RESULTS = os.path.join(BASE, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)

LENGTHS = (9, 13, 16, 20, 26, 35, 45, 55, 70)
CA_STEP_LO, CA_STEP_HI = 3.5, 4.1


def _floor(ca, phi, psi):
    """CA-RMSD between a chain and the ideal-geometry chain built from its own torsions."""
    built = geo.build_backbone(np.nan_to_num(phi), np.nan_to_num(psi))
    m = geo.kabsch_superpose(built["CA"], ca)
    return float(geo.rmsd(m, ca))


def _usable(ca, phi, psi, s, L):
    e = s + L
    step = np.linalg.norm(np.diff(ca[s:e], axis=0), axis=1)
    if ((step <= CA_STEP_LO) | (step >= CA_STEP_HI)).any():
        return False
    return bool(np.all(np.isfinite(phi[s + 1:e])) and np.all(np.isfinite(psi[s:e - 1])))


def run(n_prot=400, stride=7, seed=0, verbose=True):
    """Within-molecule length sweep over long proteins from `prots/`."""
    rng = np.random.default_rng(seed)
    paths = sorted(glob.glob(os.path.join(PROT_DIR, "*.pdb")))
    rng.shuffle(paths)
    acc = {L: [] for L in LENGTHS}
    per_prot = []
    used = 0
    t0 = time.time()
    for p in paths:
        if used >= n_prot:
            break
        try:
            seq, coords, phi, psi = geo.native_coords_from_pdb(p)
        except Exception:
            continue
        ca = np.asarray(coords["CA"], float)
        n = len(seq)
        if n < max(LENGTHS) + 2 or len(ca) != n or "X" in seq:
            continue
        phi = np.asarray(phi, float); psi = np.asarray(psi, float)
        row = {"pdb": os.path.basename(p)[:-4].upper(), "n": n}
        ok_any = False
        for L in LENGTHS:
            vals = []
            for s in range(1, n - L, stride):
                if _usable(ca, phi, psi, s, L):
                    vals.append(_floor(ca[s:s + L], phi[s:s + L], psi[s:s + L]))
            if vals:
                row[str(L)] = float(np.mean(vals))
                row[f"{L}_nwin"] = len(vals)
                acc[L].append(float(np.mean(vals)))
                ok_any = True
        if ok_any:
            per_prot.append(row)
            used += 1
            if verbose and used % 100 == 0:
                print(f"  {used}/{n_prot} proteins ({time.time()-t0:.0f}s)", flush=True)
    # PAIRED over molecules: only proteins measured at EVERY length enter the curve, so the
    # length axis is not confounded by which molecules were long enough.
    full = [r for r in per_prot if all(str(L) in r for L in LENGTHS)]
    curve = {str(L): dict(n=len(full),
                          mean=float(np.mean([r[str(L)] for r in full])),
                          median=float(np.median([r[str(L)] for r in full])),
                          se=float(np.std([r[str(L)] for r in full], ddof=1)
                                   / max(1, len(full)) ** 0.5))
             for L in LENGTHS} if full else {}
    out = dict(n_proteins=used, n_paired=len(full), lengths=list(LENGTHS),
               curve_paired=curve, rows=per_prot)
    path = os.path.join(RESULTS, "L1b_floor_vs_length.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps(dict(n_proteins=used, n_paired=len(full),
                              curve_paired=curve), indent=1))
        print("wrote", path, flush=True)
    return out


def selftest(verbose=True):
    """The estimator must return ~0 at every length on chains that ARE ideal-geometry.

    This test CAN fail: an estimator whose value grew with L for any reason other than the
    chain's own departure from ideal geometry would produce a rising curve on this input.
    """
    rng = np.random.default_rng(1)
    worst = 0.0
    curve = {}
    for L in LENGTHS:
        vals = []
        for _ in range(12):
            phi = rng.uniform(-np.pi, np.pi, L)
            psi = rng.uniform(-np.pi, np.pi, L)
            bb = geo.build_backbone(phi, psi)
            ca = np.asarray(bb["CA"], float)
            ph2, ps2 = geo.extract_torsions(bb["N"], bb["CA"], bb["C"])
            # Interior torsions are recoverable; termini are not defined from coordinates,
            # so they are taken from the generator, which is what the floor estimator does
            # for a native too.
            ph2 = np.where(np.isfinite(ph2), ph2, phi)
            ps2 = np.where(np.isfinite(ps2), ps2, psi)
            vals.append(_floor(ca, ph2, ps2))
        curve[str(L)] = float(np.mean(vals))
        worst = max(worst, float(np.max(vals)))
    if verbose:
        print(json.dumps(curve, indent=1))
        print(f"worst {worst:.2e}")
    assert worst < 1e-6, f"estimator is not exact on ideal-geometry chains: {worst}"
    # And it must be SENSITIVE: perturbing a single bond angle must move it.
    phi = rng.uniform(-np.pi, np.pi, 55); psi = rng.uniform(-np.pi, np.pi, 55)
    bb = geo.build_backbone(phi, psi)
    ca = np.asarray(bb["CA"], float).copy()
    ca[30] += 0.5
    assert _floor(ca, phi, psi) > 1e-3, "estimator is blind to a real displacement"
    print("SELFTEST PASS: exact on the manifold, sensitive off it")
    return curve


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "run"
    {"run": run, "selftest": selftest}[stage]()
