"""s24/resid0.py -- B-1.  WHAT IS ACTUALLY IN THE TORSIONAL RESIDUAL CHANNEL, AND WHAT KIND OF
ACCURACY BUYS IT.  NO MODEL, NO TRAINING, NO CORPUS -- run while Lane A's corpus does not exist.

THE QUESTION.  Lane B's architecture is `retrieved candidate + learned stochastic (dphi,dpsi) ->
new candidate`.  Before training anything it is worth knowing whether that channel HAS headroom
through the deployed readout, and -- much more importantly -- WHICH COMPONENT of a residual carries
it, because that decides what kind of model is even the right object.

THE DERIVATION UNDER TEST, written before the run.  Write a retained member as `w_k = t + e_k` and
its residual as `d_k`.  The readout is a uniform mean of 75 members, so to first order it sees
`mean_k d_k` and nothing else.  Then a residual whose per-member variation is i.i.d. contributes ~0
to the emitted cloud -- it cancels exactly as s23 L9 says the idiosyncratic 32% of member error
cancels -- and only the COMMON component moves the answer.  If that holds, a stochastic generator's
stochasticity cannot itself move the mean and the lane reduces to whether the CONDITIONAL MEAN
residual points against the 68% shared bias.  That is a claim about my own architecture and it is
tested here rather than asserted.

ARMS, n=126, point cloud everywhere, the SAME 75 members in every arm (selected once, on the real
windows, with the shipped functional, before any residual exists -- so no arm can move the answer
by changing membership).

    P0     incumbent, real windows                                        reference (3.0483)
    P0R    the same 75 rebuilt from their own torsions, no residual        THE BASELINE
    ORAC(a)  phi_k + a * d*_k      d*_k = wrap(native - member)            ORACLE
    COMM(a)  phi_k + a * dbar      dbar = circular mean of d*_k over 75    ORACLE   <- decides it
    CALT(a)  phi_k + a * wrap(native - circmean(phi_k))                    ORACLE   normalisation fork
    IDIO(a)  phi_k + a * wrap(d*_k - dbar)                                 ORACLE
    COH(s)   d*_k + eps,  ONE eps shared by all 75 members                 coherent  wrongness
    IID(s)   d*_k + eps,  eps drawn INDEPENDENTLY per member               i.i.d.    wrongness
    SHUFc/i(-)  d*_k with its RESIDUE POSITIONS permuted                   ZERO-INFORMATION CONTROL

`COH`/`IID` are the matched-accuracy corruption pair that project memory makes MANDATORY for every
learned corrector in this project (`error-coherence-decides-correctors`: at an identical 0.688 sign
accuracy, coherent mistakes emitted +0.31 A and i.i.d. mistakes -0.14 A).  They are matched on
per-member angular accuracy and differ only in whether the mistakes are shared.

`SHUF` is the plausible zero-information control IN THE OPERATOR'S OWN SPACE
(`zero-information-control-must-be-plausible`): the oracle residual's own magnitudes, permuted
across residue positions, so the per-target scale and the marginal magnitude distribution are
preserved and only the positional information is destroyed.  NOT a Gaussian of arbitrary width and
NOT "no residual", neither of which is size-matched to the arm it controls.

OPERATOR FORKS -- see `s24/PREREG_B.md` SS B-1, six axes, each naming the alternative not taken.
Enumerated by me, who has a stake in this lane; declared as a weakening, and Lane E owns the
independent enumeration.

Native torsions and native coordinates are ORACLE.  They are used for EVALUATION and for pricing a
ceiling.  Nothing here is trained and nothing here is a system result.
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

ALPHA = [0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
SIGMA = [10.0, 20.0, 30.0, 45.0, 60.0, 90.0]
OUT = os.path.join(RL.RES, "resid0.json")


def _save(o):
    t = OUT + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, OUT)


def arms_for(r, rng):
    """Every arm's (phi, psi) for one target.  Returns {name: (phi (m,n), psi (m,n))}."""
    P, S = r["PHI"], r["PSI"]
    m, n = P.shape
    dP = RL.wrap(r["nphi"][None, :] - P)          # d*_k, ORACLE
    dS = RL.wrap(r["npsi"][None, :] - S)
    bP = RL.circ_mean(dP, 0); bS = RL.circ_mean(dS, 0)            # dbar
    iP = RL.wrap(dP - bP[None, :]); iS = RL.wrap(dS - bS[None, :])  # d*_k - dbar
    # normalisation fork: the residual OF the circular mean, not the circular mean OF the residuals
    aP = RL.wrap(r["nphi"] - RL.circ_mean(P, 0)); aS = RL.wrap(r["npsi"] - RL.circ_mean(S, 0))

    out = {"P0R": (P, S)}
    for a in ALPHA:
        if a == 0.0:
            continue
        out["ORAC%.2f" % a] = (RL.wrap(P + a * dP), RL.wrap(S + a * dS))
        out["COMM%.2f" % a] = (RL.wrap(P + a * bP[None, :]), RL.wrap(S + a * bS[None, :]))
        out["CALT%.2f" % a] = (RL.wrap(P + a * aP[None, :]), RL.wrap(S + a * aS[None, :]))
        out["IDIO%.2f" % a] = (RL.wrap(P + a * iP), RL.wrap(S + a * iS))
    for s in SIGMA:
        sd = np.deg2rad(s)
        eC = rng.normal(0, sd, (1, n)); fC = rng.normal(0, sd, (1, n))
        eI = rng.normal(0, sd, (m, n)); fI = rng.normal(0, sd, (m, n))
        out["COH%02d" % s] = (RL.wrap(P + dP + eC), RL.wrap(S + dS + fC))
        out["IID%02d" % s] = (RL.wrap(P + dP + eI), RL.wrap(S + dS + fI))
    # zero-information control: the oracle residual's own values, positions permuted
    pc = rng.permutation(n)
    out["SHUFc"] = (RL.wrap(P + dP[:, pc]), RL.wrap(S + dS[:, pc]))
    pi = np.argsort(rng.random((m, n)), axis=1)
    out["SHUFi"] = (RL.wrap(P + np.take_along_axis(dP, pi, 1)),
                    RL.wrap(S + np.take_along_axis(dS, pi, 1)))
    return out


def run():
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        r = RL.retained(pdb)
        nat = r["nat"]
        rng = SD.stable_rng("s24B_resid0", pdb)
        rec = {"pdb": pdb, "n": r["n"], "fold": r["fold"]}
        # P0: the real-window incumbent, for the pinned constant only
        C0 = RL.readout(r["W"])
        rec["P0"] = float(I.ca_rmsd(C0, nat))
        e0 = RL.bias(C0, nat)
        rec["res_mag_deg"] = float(np.rad2deg(np.abs(RL.wrap(r["nphi"][None] - r["PHI"]))).mean())

        A = arms_for(r, rng)
        C0R = RL.emit(*A["P0R"])
        e0R = RL.bias(C0R, nat)                    # the UNMODIFIED rebuilt cloud's bias
        n0R = float(np.linalg.norm(e0R))
        for name, (ph, ps) in A.items():
            C = RL.emit(ph, ps)
            e = RL.bias(C, nat)
            rec[name] = float(I.ca_rmsd(C, nat))
            rec["cos_" + name] = RL.cos(e, e0)     # vs the REAL-WINDOW incumbent (the L2 spec)
            #: the coordinator's L3 measurement -- vs the arm's OWN unmodified baseline cloud.
            #: cos ~ 1 with a falling norm = the SAME error, rescaled.  cos falling = corrected.
            rec["cosR_" + name] = RL.cos(e, e0R)
            rec["nrm_" + name] = float(np.linalg.norm(e)) / n0R if n0R > 0 else float("nan")
        # validity + collapse on the two arms that matter, kept cheap
        for name in ("P0R", "ORAC0.30", "COMM0.30"):
            ph, ps = A[name]
            v = RL.validity_audit(ph, ps)
            rec["val_" + name] = {k: v[k] for k in ("rama_favoured", "clash_per_struct", "ca_ca_mean")}
            rec["col_" + name] = RL.collapse_audit(ph, ps)
        rows.append(rec)
        del r, A
        if (c + 1) % 10 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    keys = ["P0", "P0R", "SHUFc", "SHUFi"] + \
           ["%s%.2f" % (p, a) for p in ("ORAC", "COMM", "CALT", "IDIO") for a in ALPHA if a] + \
           ["%s%02d" % (p, s) for p in ("COH", "IID") for s in SIGMA]
    ok = len(rows) == len(tg) and all(all(k in r for k in keys) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "keys": keys, "alpha": ALPHA, "sigma": SIGMA, "topm": RL.TOPM})
    report(rows)


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    base = g("P0R")
    n = len(rows)

    def line(k, lab):
        s = RL.stats(g(k) - base, fold)
        print("    %-21s %7.4f %+8.4f MDE %.4f %5.2fx fold[%+.3f,%+.3f] %3dW/%3dL  %+.3f %+.3f %.3f  %s"
              % (lab, g(k).mean(), s["mean"], s["mde"], s["eff_over_mde"],
                 s["ci_fold"][0], s["ci_fold"][1], s["W"], s["L"],
                 np.nanmean(g("cos_" + k)), np.nanmean(g("cosR_" + k)), np.nanmean(g("nrm_" + k)),
                 RL.verdict(s)))

    print("\nB-1 / resid0.  n = %d.  POINT CLOUD on every side.  Baseline P0R = the same 75 windows"
          % n)
    print("rebuilt from their own torsions, no residual.  Bias vectors in the NATIVE frame,")
    print("arithmetic byte-identical to s24/qmatch.py `_bias`.  ORACLE arms labelled ORACLE.")
    print("  cosA  = bias cosine vs the REAL-WINDOW incumbent  (the L2 spec axis)")
    print("  cosR  = bias cosine vs this lane's OWN unmodified P0R cloud  (the L3 axis)")
    print("  nrm   = |bias| relative to P0R's.   cosR~1 with nrm<1 means the SAME error rescaled;")
    print("          cosR falling means the error is being CORRECTED.  Only the second is a lane.\n")
    print("    %-21s %7s %8s %27s %6s %6s %5s" % ("arm", "RMSD", "vs P0R", "", "cosA", "cosR", "nrm"))
    print("    %-21s %7.4f   (the pinned constant; P0R is the baseline, not this)"
          % ("P0  real windows", g("P0").mean()))
    print("    %-21s %7.4f   %+8.4f" % ("P0R  rebuild, no resid", base.mean(), 0.0))
    print("\n  (1) ORACLE DOSE-RESPONSE -- is there headroom in the channel at all?")
    for a in ALPHA:
        if a:
            line("ORAC%.2f" % a, "ORAC  a=%.2f  ORACLE" % a)
    print("\n  (2) THE DECOMPOSITION -- which component carries it?  THIS DECIDES THE LANE.")
    for a in ALPHA:
        if a:
            line("COMM%.2f" % a, "COMM  a=%.2f  ORACLE" % a)
    print("      -- COMM = one common correction shared by all 75 members (the conditional mean)")
    for a in ALPHA:
        if a:
            line("IDIO%.2f" % a, "IDIO  a=%.2f  ORACLE" % a)
    print("      -- IDIO = the member-specific part only, mean-free by construction")
    print("\n  (3) NORMALISATION FORK: residual-of-the-circular-mean vs circular-mean-of-residuals")
    for a in ALPHA:
        if a:
            line("CALT%.2f" % a, "CALT  a=%.2f  ORACLE" % a)
    print("\n  (4) THE ACCURACY LADDER, and the mandatory matched-accuracy coherence null.")
    print("      Both arms carry the FULL oracle direction, corrupted at the same per-member")
    print("      angular accuracy; they differ ONLY in whether the mistakes are shared.")
    for s in SIGMA:
        line("COH%02d" % s, "COH  sigma=%2d deg" % s)
    for s in SIGMA:
        line("IID%02d" % s, "IID  sigma=%2d deg" % s)
    print("\n  (5) ZERO-INFORMATION CONTROL in the operator's own space (positions permuted)")
    line("SHUFc", "SHUF coherent")
    line("SHUFi", "SHUF i.i.d.")

    print("\n  mean |oracle phi residual| = %.1f deg   (the size of what a model must predict)"
          % g("res_mag_deg").mean())
    print("\n  VALIDITY (BY CONSTRUCTION: omega dev 0, cis 0, bond len/ang dev 0 -- ideal builder)")
    for k in ("P0R", "ORAC0.30", "COMM0.30"):
        v = [r["val_" + k] for r in rows]
        print("    %-12s rama favoured %.3f   clashes/struct %.2f   Ca-Ca %.3f A"
              % (k, np.mean([x["rama_favoured"] for x in v]),
                 np.mean([x["clash_per_struct"] for x in v]),
                 np.mean([x["ca_ca_mean"] for x in v])))
    print("\n  MODE-COLLAPSE AUDIT of the 75-member sets (B=75 per target)")
    for k in ("P0R", "ORAC0.30", "COMM0.30"):
        v = [r["col_" + k] for r in rows]
        print("    %-12s unique %5.1f/75  dup %.3f  ESS %5.1f  max-mode %.3f  pairRMSD %.3f  tors-ent %.3f"
              % (k, np.mean([x["n_unique"] for x in v]), np.mean([x["dup_frac"] for x in v]),
                 np.mean([x["ess"] for x in v]), np.mean([x["mode_occ_max"] for x in v]),
                 np.mean([x["pair_rmsd_mean"] for x in v]),
                 np.mean([x["tors_entropy"] for x in v])))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
