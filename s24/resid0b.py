"""s24/resid0b.py -- B-1b.  THE CLOSING MEASUREMENT: AT ACHIEVABLE TORSION ACCURACY, IS THERE ANY
RESIDUAL AMPLITUDE THAT BEATS DOING NOTHING?

WHY THIS RUN EXISTS.  B-1 measured the accuracy ladder at FULL residual amplitude and the
pre-registered falsifier fired: carrying the whole oracle direction corrupted to sigma=30 deg is
already worse than P0R when the mistakes are shared.  That is not yet a fair closing test, because
no sensible model applies full amplitude under uncertainty -- it shrinks toward the prior, and
project memory is explicit that collapsing the posterior costs -2.253 A while being confidently
wrong costs 2-3x being absent (`torsion-restraints-reach-the-target`).  So the honest question is
the JOINT one: over the whole (amplitude, accuracy) grid, does ANY amplitude pay at the accuracy
this project can actually reach?

THE COHERENCE PARAMETER, and why it is the axis that matters here.  A residual model conditioned on
TARGET-LEVEL information -- sequence embedding, distogram, fold -- makes TARGET-LEVEL mistakes, and
those mistakes are shared by every one of the 75 members by construction.  So a learned residual's
error is COHERENT (rho -> 1) in exactly the sense project memory means, and the i.i.d. arm is the
optimistic bound that only a residual driven by per-member information could reach.  B-1b sweeps
rho so the lane is priced across that whole range instead of at its two endpoints.

    eps_k = sqrt(rho) * eps_shared + sqrt(1 - rho) * eps_k       (variance-matched at every rho)
    arm   = member torsions + alpha * (d*_k + eps_k)

ACHIEVABLE ACCURACY IS A NUMBER THIS PROJECT ALREADY OWNS, and it is the reason this run can close
the lane rather than merely describe it: sigma = 69.7 deg over all determined angles, with phi
predicted at 36.1 deg against 36.4 deg by a SEQUENCE-BLIND marginal
(`phi-carries-no-sequence-signal`).  sigma = 30 is far better than anything measured on this
corpus; sigma = 70 is where a real model sits.

LABELLING.  The best alpha at each (rho, sigma) is chosen ON THE SAME 126 TARGETS IT IS SCORED ON.
It is therefore an ORACLE amplitude and an UPPER BOUND, labelled ORACLE at every appearance, and no
number in this file is a system result.  That is deliberate: an upper bound is what closes a lane.

FORKS -- inherited unchanged from `s24/PREREG_B.md` SS B-1; this run adds no new operator choice
except the amplitude grid, whose alternative not taken is a single fixed amplitude, which would
have let a bad choice of alpha masquerade as a closed channel.
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

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
import s24.residlib as RL                 # noqa: E402

RHO = [0.0, 0.5, 1.0]
SIGMA = [30.0, 45.0, 60.0, 70.0]
ALPHA = [0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
OUT = os.path.join(RL.RES, "resid0b.json")


def _save(o):
    t = OUT + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, OUT)


def key(rho, s, a):
    return "r%.1f_s%02d_a%.2f" % (rho, s, a)


def run():
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        r = RL.retained(pdb)
        nat = r["nat"]
        rng = SD.stable_rng("s24B_resid0b", pdb)
        P, S = r["PHI"], r["PSI"]
        m, n = P.shape
        dP = RL.wrap(r["nphi"][None, :] - P)                # ORACLE direction
        dS = RL.wrap(r["npsi"][None, :] - S)
        rec = {"pdb": pdb, "n": r["n"], "fold": r["fold"]}
        C0R = RL.emit(P, S)
        rec["P0R"] = float(I.ca_rmsd(C0R, nat))
        e0R = RL.bias(C0R, nat); n0R = float(np.linalg.norm(e0R))

        for s in SIGMA:
            sd = np.deg2rad(s)
            shP = rng.normal(0, sd, (1, n)); shS = rng.normal(0, sd, (1, n))
            idP = rng.normal(0, sd, (m, n)); idS = rng.normal(0, sd, (m, n))
            for rho in RHO:
                eP = np.sqrt(rho) * shP + np.sqrt(1 - rho) * idP     # variance-matched at every rho
                eS = np.sqrt(rho) * shS + np.sqrt(1 - rho) * idS
                for a in ALPHA:
                    C = RL.emit(RL.wrap(P + a * (dP + eP)), RL.wrap(S + a * (dS + eS)))
                    k = key(rho, s, a)
                    rec[k] = float(I.ca_rmsd(C, nat))
                    e = RL.bias(C, nat)
                    rec["cosR_" + k] = RL.cos(e, e0R)
                    rec["nrm_" + k] = float(np.linalg.norm(e)) / n0R if n0R > 0 else float("nan")
        rows.append(rec)
        del r
        if (c + 1) % 10 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    keys = ["P0R"] + [key(p, s, a) for p in RHO for s in SIGMA for a in ALPHA]
    ok = len(rows) == len(tg) and all(all(k in r for k in keys) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "keys": keys,
           "rho": RHO, "sigma": SIGMA, "alpha": ALPHA, "topm": RL.TOPM})
    report(rows)


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    g = lambda k: np.array([r[k] for r in rows], float)         # noqa: E731
    base = g("P0R")
    print("\nB-1b / resid0b.  n = %d.  POINT CLOUD.  Baseline P0R = %.4f (same 75 windows rebuilt,"
          % (len(rows), base.mean()))
    print("no residual).  arm = member + alpha*(ORACLE direction + noise at accuracy sigma), with")
    print("noise coherence rho: rho=1 all 75 members share the mistake (what a TARGET-CONDITIONED")
    print("model does by construction), rho=0 independent per member.  ALL ARMS ORACLE-DIRECTED.")
    print("The alpha reported is chosen on the same 126 targets it is scored on -> ORACLE, an")
    print("UPPER BOUND.  Achievable accuracy on this corpus is sigma = 69.7 deg.\n")

    for rho in RHO:
        lab = {0.0: "rho=0.0  i.i.d. mistakes (optimistic bound)",
               0.5: "rho=0.5  half-coherent",
               1.0: "rho=1.0  SHARED mistakes  <- what a target-conditioned model does"}[rho]
        print("  %s" % lab)
        print("    %-8s %s" % ("sigma", "  ".join("a=%-5.2f" % a for a in ALPHA)))
        for s in SIGMA:
            cells = []
            for a in ALPHA:
                d = g(key(rho, s, a)) - base
                cells.append("%+7.3f" % d.mean())
            print("    %-8s %s" % ("%d deg" % s, "  ".join(cells)))
        print("    %-8s ORACLE-best alpha per sigma:" % "")
        for s in SIGMA:
            ds = [(g(key(rho, s, a)) - base, a) for a in ALPHA]
            d, a = min(ds, key=lambda x: x[0].mean())
            st = RL.stats(d, fold)
            print("      sigma=%2d  best a=%.2f  %+.4f  SE %.4f MDE %.4f  %.2fx  fold[%+.3f,%+.3f]"
                  "  %3dW/%3dL  cosR %+.3f nrm %.3f  %s  ORACLE"
                  % (s, a, st["mean"], st["se"], st["mde"], st["eff_over_mde"],
                     st["ci_fold"][0], st["ci_fold"][1], st["W"], st["L"],
                     np.nanmean(g("cosR_" + key(rho, s, a))), np.nanmean(g("nrm_" + key(rho, s, a))),
                     RL.verdict(st)))
        print()

    print("  READ THIS AGAINST sigma = 69.7 deg, the accuracy this project actually achieves,")
    print("  and against the fact that a target-conditioned residual's mistakes are rho ~ 1.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
