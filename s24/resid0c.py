"""s24/resid0c.py -- B-1c.  THE LAST DOOR: DOES AN ORACLE PER-RESIDUE CONFIDENCE GATE RESCUE THE
RESIDUAL AT ACHIEVABLE ACCURACY?

WHY.  B-1b corrupted the oracle direction with ISOTROPIC noise -- every residue equally accurate.
A real torsion model is heteroscedastic: confident and right on some residues, lost on others, and
project memory records that a confidence gate is the one construction that has ever partly rescued
a torsion channel at this length (`phi-carries-no-sequence-signal`: a PERFECT gate buys 0.5 A of a
0.94 A deficit on direct build).  Refusing to test the gated version would leave my own lane closed
on a weaker form than it deserves.

THE ARM.  At accuracy sigma with SHARED mistakes (rho = 1, what a target-conditioned model does),
apply the residual ONLY at residues where the realised error is small:

    gate_i = 1[ |eps_i| < tau ]        alpha_i = alpha * gate_i

`eps_i` is the realised corruption, so the gate is ORACLE -- it knows exactly which residues the
model got right.  No native-free confidence estimate can beat it.  tau = 180 deg is the ungated
B-1b arm and is carried as the internal control.  ORACLE at every appearance; an UPPER BOUND, and
the alpha and tau are chosen on the same 126 targets they are scored on, which loosens the bound
further in the arm's own favour.

FORKS inherited from `s24/PREREG_B.md` SS B-1.  The one addition is the gate, whose alternative not
taken is a soft confidence WEIGHT rather than a hard mask; the hard mask is the more favourable of
the two here because it cannot leak amplitude into the residues the model got wrong.
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

SIGMA = [45.0, 70.0]
TAU = [15.0, 30.0, 45.0, 180.0]
ALPHA = [0.1, 0.3, 0.5, 1.0]
OUT = os.path.join(RL.RES, "resid0c.json")


def key(s, tau, a):
    return "s%02d_t%03d_a%.2f" % (s, tau, a)


def run():
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        r = RL.retained(pdb); nat = r["nat"]
        rng = SD.stable_rng("s24B_resid0c", pdb)
        P, S = r["PHI"], r["PSI"]
        m, n = P.shape
        dP = RL.wrap(r["nphi"][None, :] - P); dS = RL.wrap(r["npsi"][None, :] - S)
        rec = {"pdb": pdb, "n": r["n"], "fold": r["fold"]}
        rec["P0R"] = float(I.ca_rmsd(RL.emit(P, S), nat))
        for s in SIGMA:
            sd = np.deg2rad(s)
            eP = rng.normal(0, sd, (1, n)); eS = rng.normal(0, sd, (1, n))   # rho = 1, shared
            for tau in TAU:
                tr = np.deg2rad(tau)
                gP = (np.abs(eP) < tr).astype(float)      # ORACLE gate
                gS = (np.abs(eS) < tr).astype(float)
                for a in ALPHA:
                    C = RL.emit(RL.wrap(P + a * gP * (dP + eP)), RL.wrap(S + a * gS * (dS + eS)))
                    rec[key(s, tau, a)] = float(I.ca_rmsd(C, nat))
                rec["cov_" + key(s, tau, 0.1)] = float((gP.mean() + gS.mean()) / 2)
        rows.append(rec); del r
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
    keys = ["P0R"] + [key(s, tau, a) for s in SIGMA for tau in TAU for a in ALPHA]
    ok = len(rows) == len(tg) and all(all(k in r for k in keys) for r in rows)
    tmp = OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "keys": keys,
                   "sigma": SIGMA, "tau": TAU, "alpha": ALPHA}, fh)
    os.replace(tmp, OUT)
    report(rows)


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    g = lambda k: np.array([r[k] for r in rows], float)        # noqa: E731
    base = g("P0R")
    print("\nB-1c / resid0c.  n = %d.  POINT CLOUD.  Baseline P0R = %.4f.  rho = 1 (SHARED mistakes,"
          % (len(rows), base.mean()))
    print("what a target-conditioned model does).  Gate is ORACLE: it keeps exactly the residues")
    print("where the realised error is under tau.  tau=180 is ungated.  ALL ARMS ORACLE.\n")
    for s in SIGMA:
        print("  sigma = %d deg" % s)
        print("    %-24s %s" % ("gate", "  ".join("a=%-5.2f" % a for a in ALPHA)))
        for tau in TAU:
            cov = g("cov_" + key(s, tau, 0.1)).mean()
            cells = ["%+7.3f" % (g(key(s, tau, a)) - base).mean() for a in ALPHA]
            print("    tau=%3d deg (keeps %4.1f%%)  %s" % (tau, 100 * cov, "  ".join(cells)))
        best = min([(g(key(s, tau, a)) - base, tau, a) for tau in TAU for a in ALPHA],
                   key=lambda x: x[0].mean())
        d, tau, a = best
        st = RL.stats(d, fold)
        print("    ORACLE-best cell: tau=%d a=%.2f  %+.4f  SE %.4f MDE %.4f  %.2fx"
              "  fold[%+.3f,%+.3f]  %3dW/%3dL  %s  ORACLE"
              % (tau, a, st["mean"], st["se"], st["mde"], st["eff_over_mde"],
                 st["ci_fold"][0], st["ci_fold"][1], st["W"], st["L"], RL.verdict(st)))
        print()


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
