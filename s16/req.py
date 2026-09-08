"""s16/req.py -- THE REQUIREMENT CURVE.  How good must a direction be before steering pays?

Motivated and written BEFORE the n = 126 flagship numbers returned (the run was in flight;
see `s16/PREREG_steer.md` and ledger L9).  It exists so that a null in `s16/steer.py` is not
merely a null: whatever the native-free direction does, this module says what a direction
would have to BE for the mechanism to work at all, and whether the loud projection helps a
direction that is already good.

CONSTRUCTION.  For each target we interpolate between the native-free direction `u` and the
ORACLE error `e`, at unit norm:

    v(a) = normalise( a * e_hat + (1 - a) * u_hat ),      a in [0, 1]

and rescale to ||e||, so every rung is a step of the SAME LENGTH that differs ONLY in its
cosine with the truth.  That isolates DIRECTION QUALITY from step length, which is the whole
question.  The realised cosine is measured, never assumed -- the map from `a` to cos is
non-linear and target-dependent, exactly as the mixing weight is in `s16/qphase.py`.

EVERY ARM HERE IS AN ORACLE DIAGNOSTIC.  `v(a)` reads the native for any a > 0.  This is a
REQUIREMENT CURVE and a CEILING, in the same class as Sprint 15's "perfect distance
knowledge is worth 0.611 A".  No number in this module may be quoted as a predictive result,
and the a = 0 rung -- which is native-free -- is the only rung that coincides with anything
achievable today.

WHAT IT ANSWERS, in one table:

  1. the cosine at which a full-space step first beats the plain fit out of fold;
  2. whether the LOUD PROJECTION raises or lowers that threshold -- the flagship's actual
     hypothesis, asked of directions good enough for the answer to be visible;
  3. how far the curve is from the |cos| ~ 0.36 that the native-free channel supplies.

If the loud projection lowers the required cosine, the mechanism is real and the native-free
channel is simply too weak, which is a different and more actionable conclusion than "the
projection does not work".  If the projection RAISES the threshold at every cosine, the
projection is harmful in itself and no improvement to `u` would rescue it.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s16.steer import (pool_torsions, rmsd_of, subspaces,   # noqa: E402
                       _native_torsions)

ALPHAS = (0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0)
STEPS = (0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0)
RANKS = (4, 10)


def _unit(v):
    return v / max(float(np.linalg.norm(v)), 1e-12)


def run(targets=None, n_start=4):
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    data = C.gather(tg)

    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    rows = []
    for c, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        i, j, n, fold = d["i"], d["j"], d["n"], d["fold"]
        sd = d["sd"]; nat = d["nat"]
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)

        starts = D.starts(p, d["seq"], n, fold, n_start, SD.stable_rng(p, "s16steer"))
        best = None
        for phi0, psi0, _tag in starts:
            ph, ps, f = A.fit(dhat, sd, i, j, phi0, psi0)
            if best is None or f < best[0]:
                best = (f, ph, ps)
        _f, phi, psi = best
        th = np.concatenate([phi, psi])

        _J, V, _sv, _CA = subspaces(phi, psi)
        nphi, npsi = _native_torsions(p, nat, d["seq"], fold)
        e = A.wrap(np.concatenate([nphi, npsi]) - th)              # ORACLE
        ppl, psl = pool_torsions(p, n, fold)
        u = A.wrap(np.concatenate([ppl, psl]) - th)                # native-free
        ne = float(np.linalg.norm(e))
        eh, uh = _unit(e), _unit(u)

        base = rmsd_of(phi, psi, nat)
        cell = {}
        for a in ALPHAS:
            v = _unit(a * eh + (1.0 - a) * uh) * ne
            #: the REALISED cosine, measured rather than assumed.
            cos = float(abs(_unit(v) @ eh))
            variants = {"raw": v}
            for r in RANKS:
                Vr = V[:, :r]
                variants[f"loud{r}"] = (Vr @ Vr.T) @ v
            lad = {}
            for vn, vec in variants.items():
                lad[vn] = {}
                for s in STEPS:
                    x = th + s * vec
                    lad[vn][str(s)] = rmsd_of(x[:n], x[n:], nat)
            cell[str(a)] = {"cos": cos, "lad": lad}

        rows.append({"pdb": p, "n": int(n), "fold": int(fold), "base": base,
                     "err_norm": ne, "cos_u_e": float(abs(uh @ eh)), "cells": cell})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)
            json.dump({"n": len(rows), "rows": rows},
                      open(os.path.join(RESULTS, "req.json"), "w"))

    json.dump({"n": len(rows), "rows": rows},
              open(os.path.join(RESULTS, "req.json"), "w"))
    report(rows, folds)
    return rows


def _lfo(rows, folds, a, vn):
    """leave-fold-out step choice, the same rule the flagship ladder uses."""
    picked, vals = {}, np.zeros(len(rows))
    for f in sorted(set(folds)):
        tr = folds != f
        picked[int(f)] = min(STEPS, key=lambda s: np.mean(
            [rows[q]["cells"][str(a)]["lad"][vn][str(s)]
             for q in range(len(rows)) if tr[q]]))
    for q, r in enumerate(rows):
        vals[q] = r["cells"][str(a)]["lad"][vn][str(picked[int(folds[q])])]
    return vals, picked


def report(rows, folds):
    folds = np.asarray(folds, int)
    base = np.array([r["base"] for r in rows])
    print(f"\nn = {len(rows)}   plain fit {base.mean():.3f} A   "
          f"native-free |cos(u,e)| = {np.mean([r['cos_u_e'] for r in rows]):.3f}")
    print("\nALL ARMS BELOW ARE ORACLE DIAGNOSTICS -- the direction is built from the native.")
    print("Every rung is a step of the SAME LENGTH ||e||; only its cosine with the truth changes.\n")
    hdr = f"  {'alpha':>6}{'realised cos':>14}"
    for vn in ["raw"] + [f"loud{r}" for r in RANKS]:
        hdr += f"{vn + ' RMSD':>14}{'vs fit':>10}"
    print(hdr)
    for a in ALPHAS:
        cos = np.mean([r["cells"][str(a)]["cos"] for r in rows])
        line = f"  {a:>6.2f}{cos:>14.3f}"
        for vn in ["raw"] + [f"loud{r}" for r in RANKS]:
            v, _pk = _lfo(rows, folds, a, vn)
            line += f"{v.mean():>14.3f}{v.mean() - base.mean():>+10.3f}"
        print(line)
    print("\n  READ: the row where 'vs fit' first turns negative is the cosine a direction must"
          "\n  reach before a step of the true error's length pays. If loud4/loud10 turn negative"
          "\n  at a LOWER cosine than raw, the loud projection is doing real work and the"
          "\n  native-free channel is merely too weak; if they turn negative LATER or never, the"
          "\n  projection is harmful in itself and improving `u` would not rescue it.")


if __name__ == "__main__":
    run()
