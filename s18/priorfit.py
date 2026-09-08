"""s18/priorfit.py -- THE ARM MATH'S MECHANISM OPENS.

WHY THIS EXISTS.  The degree-1 branch is closed (F2/F4/F5), but MATH's closure named a
mechanism that points somewhere new.  On the 19-target enumerated instrument the objective
`hamil` is **75% an exactly-additive torsion prior**, and the residue-additive projection of
it correlates with the pure retrieval-pool torsion prior at **rho = 0.986** with matching
argmins.  So the "degree-1 object" whose argmin looked good was, to three nines, THE TORSION
PRIOR.

    And `s17/refine.py`'s continuum objective contains NO PRIOR TERM AT ALL.

That is why nothing transferred: the truncation had nothing to reduce to.  It also says the
continuum objective is missing the exact term that was carrying the enumerated instrument's
argmin quality.  This module adds it back.

    THE QUESTION: does a native-free per-residue torsion prior, added to the deployed
    distance objective, rescue refinement from the +0.561 A it currently emits?

THE MOTIVATING DEFECT.  Sprint 17 measured the unregularised fit to be a LARGE-DISPLACEMENT
operation -- 87% of the displacement of a random torsion assignment.  A prior term penalises
exactly that.  This is not "optimise the current objective harder" (a standing refutation);
it is a different objective, and the term added is the best torsion channel measured anywhere
in this project (`s14/retprior`: phi 33.6 deg, psi 59.2 deg, beating the trained leave-fold-out
sequence predictor).

THE OBJECTIVE.  Both terms standardised by their own count so `lam` is interpretable:

    E(theta) = (1/np) * sum_p ((d_p - dhat_p)/sd_p)^2
             + lam * (1/nr) * sum_i [ kphi_i (1 - cos(phi_i - mphi_i))
                                    + kpsi_i (1 - cos(psi_i - mpsi_i)) ]

The prior is a per-residue von Mises fitted to the RETRIEVED WINDOWS: circular mean for the
location, and concentration from the mean resultant length via the standard approximation
kappa ~ R(2 - R^2)/(1 - R^2).  Everything is native-free; the native is read only to score.

ARMS.  `lam = 0` must reproduce `s17/refine.py::refine_full` -- a built-in reproduction check.

    avg           the coordinate average of the shipped top-75          the start, ~3.048 A
    proj          its ideal-geometry projection                         the torsion start
    prior{lam}    distance + lam * retrieval prior                      the hypothesis
    prioronly     the prior alone, no distance term                     what the prior knows
    CONTROLS, both at every lam
    helix{lam}    IDENTICAL functional form, mu = ideal alpha-helix      ZERO-INFORMATION prior
    shuf{lam}     the per-residue (mu, kappa) tuples PERMUTED across i   matched-marginal

THE CONTROL THAT DECIDES IT, and the reason it is not optional.  The programme has already
recorded that on this instrument **a zero-information constant alpha-helix beats the random
control**, so "beats random" proves nothing here.  `helix` carries the same functional form,
the same stiffness distribution and the same regularising pull, and knows NOTHING about which
residue is which.  If `helix` captures the effect, the gain is REGULARISATION -- shrinking a
large-displacement fit -- and not the prior's positional information.  `shuf` separates the
same two things a second way.  Both must lose for the prior's information to be credited.

PRE-REGISTRATION, written before the run.

  HYPOTHESIS.  Some lam > 0 beats the coordinate average at n = 126.

  EXPECTED, honestly.  I expect a real effect and I expect `helix` to take most of it.  The
  s17 diagnosis is a displacement pathology, and displacement is exactly what a zero-information
  prior fixes.  The live question is whether ANY of the gain survives the helix control.

  SUCCESS.  Some `prior{lam}` beats `avg` at n = 126, target as the unit, paired fold-clustered
  interval excluding zero, AND beats `helix{lam}` and `shuf{lam}` at the SAME lam with intervals
  excluding zero.

  FALSIFIER.  If the best `prior{lam}` is within the instrument's 0.084 A MDE of `helix{lam}`,
  the prior's positional information contributes nothing and the finding is "the fit needed
  regularising", which is a smaller and different claim that must be reported as such.

  THE TRAP I AM WATCHING FOR.  `lam` -> infinity walks toward `prioronly`, which is a pool
  summary and not a fit at all.  A win at the largest lam is a statement that the distance term
  should be DOWN-WEIGHTED, not that the prior helps -- read the lam curve, never one rung.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from core import project as pj              # noqa: E402
from s12 import instrument as I             # noqa: E402
from s14 import avgspace as AV              # noqa: E402
from s14 import ladder as LD                # noqa: E402
from s14 import retprior as RP              # noqa: E402
from s15 import align_lib as A              # noqa: E402
from s15 import distcal as C                # noqa: E402
from s15 import seed as SD                  # noqa: E402

LAMS = (0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0)


def von_mises(pdb, arm="top75"):
    """Native-free per-residue (mu, kappa) for phi and psi from the RETRIEVED windows.

    kappa from the mean resultant length by the standard approximation; clipped so a residue
    the pool disagrees about cannot contribute an unbounded restraint.
    """
    PHI, PSI, sim = RP.windows(pdb, arm)
    out = []
    for X in (PHI, PSI):
        mu = RP.circ_mean(X, w=sim, axis=0)
        R = np.clip(RP.circ_R(X, w=sim, axis=0), 1e-6, 0.999)
        k = R * (2.0 - R * R) / (1.0 - R * R)
        out.append((np.asarray(mu, float), np.clip(np.asarray(k, float), 0.0, 50.0)))
    return out[0], out[1]


def fit_prior(dhat, sd, i, j, phi0, psi0, mphi, kphi, mpsi, kpsi, lam, maxiter=400):
    """L-BFGS on the deployed distance objective PLUS lam * the von Mises prior.

    The distance half is identical in form to `s15.align_lib.fit`'s squared branch, so
    `lam = 0` reproduces `refine_full` rather than merely resembling it.
    """
    from scipy.optimize import minimize

    n = len(phi0)
    inv = 1.0 / np.asarray(sd, float) if len(dhat) else np.zeros(0)
    #: The distance half is left UNNORMALISED so it is bit-identical to `align_lib.fit`'s
    #: `kind="squared"` branch (rho = r^2, drho = 2r) and `lam = 0` REPRODUCES `refine_full`.
    #: A per-pair rescaling was tried first and failed its own reproduction gate on 2 of 3
    #: smoke targets: dividing f by npair shrinks the gradient ~100x, `gtol` is absolute, and
    #: the non-convex fit stopped in a different basin.  The scale ratio moves onto the prior
    #: instead, so `lam` keeps its meaning -- per-pair against per-residue contribution.
    scale = float(max(len(dhat), 1)) / float(n)

    def fg(x):
        phi, psi = x[:n], x[n:]
        f = 0.0
        g = np.zeros(2 * n)
        if len(dhat):
            G = pj.frames(phi[None], psi[None])[0]
            CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
            rv = CA[i] - CA[j]
            d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
            r = (d - dhat) * inv
            f = float((r * r).sum())
            coef = ((2.0 * r * inv) / d)[:, None] * rv
            gCA = np.zeros_like(CA)
            np.add.at(gCA, i, coef)
            np.add.at(gCA, j, -coef)
            g = pj._torsion_grad(G, CA, gCA)
        if lam > 0.0:
            w = lam * scale
            dp, ds = phi - mphi, psi - mpsi
            f += w * float((kphi * (1.0 - np.cos(dp))).sum()
                           + (kpsi * (1.0 - np.cos(ds))).sum())
            g = g + np.concatenate([w * kphi * np.sin(dp), w * kpsi * np.sin(ds)])
        return f, g

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    res = minimize(fg, x0, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return res.x[:n], res.x[n:], float(res.fun)


def run(targets=None):
    #: The leave-fold-out distogram debias is ALWAYS fitted on the full target list, even when
    #: `targets` is a subset.  Fitting it on the subset silently changes `dhat` and made a
    #: 3-target smoke look like a failed reproduction gate against s17's n = 126 numbers.
    full = I.targets()
    tg = targets if targets is not None else full
    rows, t0 = [], time.time()
    data = C.gather(full)
    pdbs = [t["pdb"] for t in full]
    deb = {}
    for f in sorted({int(t["fold"]) for t in full}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
        rng = SD.stable_rng(pdb, "s18priorfit")

        W = np.asarray(AV.top75_windows(pdb)[0], float)
        P = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

        (mphi, kphi), (mpsi, kpsi) = von_mises(pdb)
        #: CONTROL -- zero-information prior: ideal alpha-helix everywhere, SAME stiffness.
        hphi = np.full(n, LD.HELIX_PHI)
        hpsi = np.full(n, LD.HELIX_PSI)
        #: CONTROL -- matched-marginal: the per-residue (mu, kappa) tuples permuted across i.
        pm = rng.permutation(n)

        e = {"avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
             "proj": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)),
             "kphi": float(kphi.mean()), "kpsi": float(kpsi.mean())}

        for lam in LAMS:
            for tag, aa in (("prior", (mphi, kphi, mpsi, kpsi)),
                            ("helix", (hphi, kphi, hpsi, kpsi)),
                            ("shuf", (mphi[pm], kphi[pm], mpsi[pm], kpsi[pm]))):
                if lam == 0.0 and tag != "prior":
                    continue                      #: lam = 0 is the same fit in all three arms
                p_, q_, f_ = fit_prior(dhat, sd, i, j, phi0, psi0,
                                       aa[0], aa[1], aa[2], aa[3], lam)
                key = f"{tag}{lam:g}"
                e[key] = float(I.ca_rmsd(I.build_ca(p_, q_), nat))
                e["obj_" + key] = f_
                e["disp_" + key] = float(np.linalg.norm(
                    A.wrap(np.concatenate([p_, q_]) - np.concatenate([phi0, psi0]))))

        #: the prior ALONE -- what the torsion channel knows with no distance term at all
        z = np.zeros(0, int)
        p_, q_, _f = fit_prior(np.zeros(0), np.ones(0), z, z, phi0, psi0,
                               mphi, kphi, mpsi, kpsi, 1.0)
        e["prioronly"] = float(I.ca_rmsd(I.build_ca(p_, q_), nat))

        rows.append({"pdb": pdb, "n": n, "fold": fold, **e})
        if (c + 1) % 5 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": False},
                      open(os.path.join(RESULTS, "priorfit.json"), "w"))
    json.dump({"rows": rows, "complete": len(rows) == len(tg)},
              open(os.path.join(RESULTS, "priorfit.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "priorfit.json")))["rows"]
    rng = SD.stable_rng("priorfit", "report")
    g = lambda k: np.array([r[k] for r in rows])            # noqa: E731
    avg = g("avg")
    print(f"\nn = {len(rows)}.  Start = the coordinate average of the shipped top-75 "
          f"({avg.mean():.3f} A).  All arms NATIVE-FREE.\n")
    print(f"  {'arm':<14}{'RMSD':>8}{'median':>9}{'vs avg [95% CI]':>26}{'W/L':>9}"
          f"{'disp(rad)':>11}")
    print(f"  {'avg (start)':<14}{avg.mean():>8.3f}{np.median(avg):>9.3f}{'--':>26}"
          f"{'--':>9}{'--':>11}")
    keys = ["proj"]
    for lam in LAMS:
        for tag in ("prior", "helix", "shuf"):
            k = f"{tag}{lam:g}"
            if k in rows[0]:
                keys.append(k)
    keys.append("prioronly")
    for k in keys:
        v = g(k)
        mu, lo, hi = _boot(v - avg, rng)
        w = int((v < avg).sum()); l = int((v > avg).sum())
        dk = "disp_" + k
        dd = f"{g(dk).mean():>11.3f}" if dk in rows[0] else f"{'--':>11}"
        tag = "  CONTROL" if k.startswith(("helix", "shuf")) else ""
        print(f"  {k:<14}{v.mean():>8.3f}{np.median(v):>9.3f}   "
              f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l}{dd}{tag}")

    print("\n  PRIOR vs ITS CONTROLS at the same lam -- this is the row that decides it:")
    for lam in LAMS:
        if lam == 0.0:
            continue
        p = g(f"prior{lam:g}")
        for tag in ("helix", "shuf"):
            k = f"{tag}{lam:g}"
            if k not in rows[0]:
                continue
            mu, lo, hi = _boot(p - g(k), rng)
            print(f"    lam={lam:<6g} prior - {tag:<6}{mu:+.3f} [{lo:+.3f},{hi:+.3f}]"
                  f"{int((p < g(k)).sum()):>6}/{int((p > g(k)).sum())}")
    print("\nREAD.  `lam=0` must reproduce s17's refine_full (3.610).  A win over `avg` that is")
    print("NOT also a win over `helix` at the same lam is REGULARISATION, not the prior's")
    print("positional information -- report it as the smaller claim.  Read the lam CURVE: a win")
    print("only at the largest lam says the DISTANCE term should be down-weighted.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
