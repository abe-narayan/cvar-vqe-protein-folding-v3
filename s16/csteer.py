"""s16/csteer.py -- COORDINATE-SPACE STEERING.  The evidence-driven pivot from the flagship.

WHY THIS EXISTS.  `s16/steer.py` refuted torsion-space steering at n = 126 and, more usefully,
said WHY: stepping the EXACT native error in torsion space realises 5.0% of the achievable
RMSD gain at half-way and 29.2% at seven-tenths.  RMSD is strongly convex -- over the first
third, non-monotonic -- in torsion displacement toward the truth.  Steering requires that
partial motion toward the truth pay partially, and in torsion space it does not.

In COORDINATE space it does, and it does so EXACTLY, by algebra rather than by hope.  Let `r`
be the residual (start minus native, after superposition) and `v` a unit direction with
`cos = v . r_hat`.  Stepping `t` along `v`:

    ||r - t v||^2 = ||r||^2 - 2 t (v . r) + t^2      minimised at   t* = v . r = c ||r||
    giving                                    ||r_new|| = ||r|| sqrt(1 - c^2)

So in coordinate space ANY direction with c > 0 pays, the optimal step is known in closed
form, and the entire question collapses to ONE measurable native-free quantity: **c**.

That also turns the sprint's headline goal into an arithmetic requirement, which is the most
useful thing this module produces:

    3.204 A -> 2.5 A  needs  c = 0.625        3.204 A -> 2.0 A  needs  c = 0.781

Sprint 15's parameter-free FUSION LAW, `d_avg ~ sqrt(r^2 - (s/2)^2)`, is the `t = ||d||/2`
special case of the same algebra, already validated to 0.162 A on 126 targets.  This module
generalises it from the fixed half-step to the optimal step, and measures `c` for every
native-free channel the programme has.

DESIGN.  Three START structures and, for each, every native-free direction available:

    starts       fit    the torsion-space distance-geometry fit          (3.647 A)
                 avg    the top-75 coordinate average, pre-projection
                 proj   that average put through the ideal-geometry projection
    directions   toward each OTHER start, toward the multi-start ensemble mean, and
                 toward the pool medoid -- all native-free, all superposed onto the start
    controls     rand   matched-norm random coordinate direction (3 draws, stable_rng)
                 ORACLE_true  the exact residual, whose optimal step is the ceiling

Step chosen leave-fold-out on the same ladder for every arm, with `s = 0` present so a
training fold can always decline to move.  The ORACLE optimal step `t* = c ||r||` is reported
beside it as a labelled ceiling, never as a result.

WHAT WOULD FALSIFY THE PIVOT.  If the measured `c` for every native-free channel is at or
below the torsion-space 0.318, coordinate-space steering buys at most
`r (1 - sqrt(1 - 0.318^2)) = 0.052 r` ~ 0.19 A from the fit and cannot reach the incumbent.
That is a real possible outcome and it is written down here before the run.

VALIDITY.  A steered coordinate set is not a valid peptide.  The `proj` readout
(`python -m s16.csteer proj`) puts the winning arm back through the ideal-geometry projection
and re-scores it, so no accuracy number is ever quoted off an unbuildable structure.
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

STEPS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0)
N_RAND = 3


def _sup(X, T):
    """superpose X onto T (both (n, 3)), returning X in T's frame."""
    return I.superpose_batch(np.asarray(X, float)[None], np.asarray(T, float))[0]


def _cos(d, r):
    nd = float(np.linalg.norm(d)); nr = float(np.linalg.norm(r))
    if nd < 1e-12 or nr < 1e-12:
        return 0.0
    return float((d.ravel() @ r.ravel()) / (nd * nr))


def run(targets=None, n_start=4):
    from s14 import avgspace as AV
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
    for c_i, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        i, j, n, fold = d["i"], d["j"], d["n"], d["fold"]
        sd = d["sd"]; nat = np.asarray(d["nat"], float); seq = d["seq"]
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)

        # ---- channel 1: the torsion-space fit, and its multi-start ensemble
        starts_t = D.starts(p, seq, n, fold, n_start, SD.stable_rng(p, "s16steer"))
        best, ens = None, []
        for phi0, psi0, _tag in starts_t:
            ph, ps, f = A.fit(dhat, sd, i, j, phi0, psi0)
            ens.append(I.build_ca(ph, ps))
            if best is None or f < best[0]:
                best = (f, ph, ps)
        CA_fit = I.build_ca(best[1], best[2])

        # ---- channel 2: the retrieval pool, its average, medoid and projection
        W = np.asarray(AV.top75_windows(p)[0], float)
        P = I.pairwise_rmsd(W)
        avg_c, b_med = I.coordinate_average(W, P)
        avg = np.asarray(avg_c, float)
        med = np.asarray(W[b_med], float)
        proj = np.asarray(I.project(avg, seq, fold)["ca"], float)
        #: the ensemble mean is taken in COORDINATE space, superposed member by member --
        #: the flagship's whole lesson is that the torsion-space mean is the wrong object.
        ens_m = np.mean([_sup(E, CA_fit) for E in ens], axis=0)

        pool = {"fit": CA_fit, "avg": avg, "proj": proj, "med": med, "ens": ens_m}
        rng = SD.stable_rng(p, "s16csteer")

        cell = {}
        for sname in ("fit", "avg", "proj"):
            S = pool[sname]
            Nsup = _sup(nat, S)                      # ORACLE: native in the start's frame
            r = Nsup - S                             # ORACLE residual
            nr = float(np.linalg.norm(r))
            base = float(I.ca_rmsd(S, nat))

            dirs = {}
            for tname, T in pool.items():
                if tname == sname:
                    continue
                dirs[tname] = _sup(T, S) - S
            for k in range(N_RAND):
                z = rng.standard_normal(S.shape)
                z -= z.mean(0)                       # a pure translation is invisible to RMSD
                nz = float(np.linalg.norm(z))
                ref = float(np.linalg.norm(dirs["avg" if sname != "avg" else "fit"]))
                dirs[f"rand{k}"] = z * (ref / max(nz, 1e-12))
            dirs["ORACLE_true"] = r

            arms = {}
            for dname, dv in dirs.items():
                nd = float(np.linalg.norm(dv))
                cc = _cos(dv, r)
                lad = {str(s): float(I.ca_rmsd(S + s * dv, nat)) for s in STEPS}
                #: the closed-form optimum, in the direction's OWN units.  ORACLE: it uses r.
                t_star = (cc * nr / nd) if nd > 1e-12 else 0.0
                arms[dname] = {
                    "cos": cc, "norm": nd, "lad": lad,
                    "t_star": float(t_star),
                    "rmsd_at_t_star": float(I.ca_rmsd(S + t_star * dv, nat)),
                    #: the algebraic prediction ||r|| sqrt(1 - c^2), in RMSD units.
                    "law_pred": float(base * np.sqrt(max(0.0, 1.0 - cc * cc))),
                }
            cell[sname] = {"base": base, "resid_norm": nr, "arms": arms}

        rows.append({"pdb": p, "n": int(n), "fold": int(fold), "cells": cell})
        if (c_i + 1) % 10 == 0:
            print(f"  {c_i+1}/{len(tg)} checkpointed", flush=True)
            json.dump({"n": len(rows), "rows": rows},
                      open(os.path.join(RESULTS, "csteer.json"), "w"))

    json.dump({"n": len(rows), "rows": rows},
              open(os.path.join(RESULTS, "csteer.json"), "w"))
    report(rows, folds)
    return rows


def _boot(diff, rng, B=4000):
    diff = np.asarray(diff, float); k = len(diff)
    m = diff[rng.integers(0, k, size=(B, k))].mean(1)
    return float(diff.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def _lfo(rows, folds, sname, dname):
    picked, vals = {}, np.zeros(len(rows))
    for f in sorted(set(folds)):
        tr = folds != f
        picked[int(f)] = min(STEPS, key=lambda s: np.mean(
            [rows[q]["cells"][sname]["arms"][dname]["lad"][str(s)]
             for q in range(len(rows)) if tr[q]]))
    for q, r in enumerate(rows):
        vals[q] = r["cells"][sname]["arms"][dname]["lad"][str(picked[int(folds[q])])]
    return vals, picked


def report(rows, folds):
    folds = np.asarray(folds, int)
    rng = SD.stable_rng("csteer", "report")
    print(f"\nn = {len(rows)}")
    print("\nIn coordinate space the payoff is exact algebra: a direction with cosine c, stepped")
    print("optimally, takes RMSD r to r*sqrt(1-c^2).  So c IS the result.  Requirement, from the")
    print("incumbent's 3.204 A:   2.5 A needs c = 0.625      2.0 A needs c = 0.781\n")

    for sname in ("fit", "avg", "proj"):
        base = np.array([r["cells"][sname]["base"] for r in rows])
        names = [k for k in rows[0]["cells"][sname]["arms"] if not k.startswith("rand")]
        names += ["rand0"]
        print(f"START = {sname}   mean {base.mean():.3f} A")
        print(f"  {'direction':<14}{'cos [95% CI]':>26}{'law r*sqrt(1-c^2)':>19}"
              f"{'ORACLE t*':>11}{'LFO step':>26}")
        for dname in names:
            cs = np.array([r["cells"][sname]["arms"][dname]["cos"] for r in rows])
            mu, lo, hi = _boot(cs, rng)
            law = np.mean([r["cells"][sname]["arms"][dname]["law_pred"] for r in rows])
            att = np.mean([r["cells"][sname]["arms"][dname]["rmsd_at_t_star"] for r in rows])
            v, pk = _lfo(rows, folds, sname, dname)
            dm, dlo, dhi = _boot(v - base, rng)
            tag = "  ORACLE" if dname.startswith("ORACLE") else ""
            print(f"  {dname:<14}{mu:>10.3f} [{lo:+.3f},{hi:+.3f}]{law:>19.3f}{att:>11.3f}"
                  f"   {dm:+.3f} [{dlo:+.3f},{dhi:+.3f}]{tag}")
        print()

    print("READ.  `law` and `ORACLE t*` are CEILINGS -- both use the native residual to place the")
    print("step.  The only deployable column is `LFO step`, whose step is trained leave-fold-out")
    print("on labels.  A gap between `law` and `LFO step` is the price of not knowing t*.")


if __name__ == "__main__":
    run()
