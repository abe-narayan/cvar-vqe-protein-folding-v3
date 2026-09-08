"""SPRINT 15, coordinator -- ROBUST restraint losses, because the violations are CONCENTRATED.

WHY THIS EXPERIMENT EXISTS, AND WHY IT SHOULD HAVE BEEN THE FIRST ONE.

`s15/feasible.py` measured, on all 126 targets and 8,549 pairs, how far the NATIVE structure sits
from its own predicted restraints, in units of the predictor's own reported uncertainty:

    median z = 1.336        max z = 7.237
    eps 0.5 -> 24.8% of pairs admit the native        eps 2.0 -> 71.2%
    eps 1.0 -> 45.8%                                  eps 3.0 -> 84.0%

Read that shape carefully, because it is not the shape of a badly calibrated predictor. A typical
restraint is only MILDLY inconsistent with the truth -- half of them are inside 1.34 sigma. The
damage is done by a MINORITY of badly wrong pairs, the worst of them off by more than SEVEN of the
predictor's own standard deviations. That is an **outlier** problem, not a calibration problem.

And squared error is the worst possible loss for an outlier problem. A pair at z = 7 contributes
**49 times** the pull of a pair at z = 1. Under weighted least squares the handful of restraints
that are catastrophically wrong dominate the fit, drag the structure away from the many restraints
that are approximately right, and do so with more force the more confident the predictor was.
K2 explains why recalibration cannot rescue this: the z-sd is nearly CONSTANT across separation
(2.19, 2.98, 2.70, 2.70, 2.67), so the `sd` is wrong by a near-uniform factor, and a uniform
rescaling of all weights leaves a weighted least-squares argmin **exactly unchanged**. It cannot
separate a mildly-wrong restraint from a catastrophically-wrong one. Only the LOSS can.

This is also not an exotic idea. Restrained refinement in NMR structure determination has handled
exactly this failure mode for thirty years, with soft-square and robust potentials and with
violation-aware reweighting, precisely because experimentally derived distance restraints contain
misassignments. This project has been using plain squared error the whole time.

THE ARMS. Every arm is the identical fit -- same starts, same analytic O(n) torsion gradient, same
objective-only multi-start selection -- differing ONLY in the loss applied to the standardised
residual r = (d_ij - dhat_ij) / sd_ij.

    squared      rho(r) = r^2                            the incumbent objective, the control
    huber        r^2 for |r| <= delta, else linear       classical bounded-influence
    soft_l1      2 (sqrt(1 + r^2) - 1)                   smooth, mild
    cauchy       log(1 + r^2)                            strong redescending influence
    welsch       1 - exp(-r^2)                           fully redescending: a bad pair is IGNORED
    trimmed      squared error over the best (1-q) of    explicit violation rejection, the
                 pairs, refitted                          direct NMR analogue

`cauchy` and `welsch` are redescending: their influence function returns toward zero for large
residuals, so a restraint that is wildly inconsistent stops pulling at all. That is the behaviour
the measured z-profile calls for, and it is also the risky one -- a redescending loss can happily
ignore a restraint that was RIGHT and the structure wrong, which is how it fails.

WHAT WOULD FALSIFY THE IDEA. If no robust loss beats `squared` on the full instrument, then the
outlier reading of K7 is wrong: the restraint errors would be diffuse rather than concentrated, and
the deficit is genuine information loss that no loss function can recover. That is a clean negative
and it is worth having, because it would close a direction that looks obviously promising.

HYPERPARAMETERS ARE CHOSEN LEAVE-FOLD-OUT. `delta`, the Cauchy/Welsch scale, and the trim fraction
`q` are selected on the four training folds by the OBJECTIVE-selected fit's restraint residual --
never by RMSD, and never on the fold being reported. An arm tuned on the targets it is reported on
is worthless and the brief forbids it. `ORACLE_scale` is included as the in-fold ceiling: it says
whether a disappointing result is a transfer failure or a functional-form failure, the distinction
Sprint 14 showed decides the interpretation of every negative.

Run:
    python -m s15.robust
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from core import project as pj               # noqa: E402
from s12 import instrument as I              # noqa: E402
from s15 import seed as SD             # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

LOSSES = ("squared", "huber", "soft_l1", "cauchy", "welsch")
SCALES = (0.5, 1.0, 2.0, 3.0)          # the loss scale, in units of the predictor's own sd
TRIMS = (0.05, 0.10, 0.20)


# ------------------------------------------------------------------------------ the losses
def rho(r, kind, s):
    """Loss and its derivative, both as functions of the standardised residual `r`.

    Returned as (value, d rho / d r). Every loss is scaled so that it agrees with `squared`
    in the small-residual limit, which keeps the arms comparable at the same `s`.
    """
    if kind == "squared":
        return r * r, 2.0 * r
    x = r / s
    if kind == "huber":
        a = np.abs(x)
        small = a <= 1.0
        v = np.where(small, x * x, 2.0 * a - 1.0)
        g = np.where(small, 2.0 * x, 2.0 * np.sign(x))
        return (s * s) * v, s * g
    if kind == "soft_l1":
        t = np.sqrt(1.0 + x * x)
        return (s * s) * 2.0 * (t - 1.0), s * (2.0 * x / t)
    if kind == "cauchy":
        return (s * s) * np.log1p(x * x), s * (2.0 * x / (1.0 + x * x))
    if kind == "welsch":
        e = np.exp(-x * x)
        return (s * s) * (1.0 - e), s * (2.0 * x * e)
    raise ValueError(kind)


def fit_robust(dhat, sd, i, j, phi0, psi0, kind="squared", s=1.0, mask=None, maxiter=400):
    """Fit torsions to restraints under a robust loss, with the exact analytic gradient.

    The chain rule is identical to plain least squares -- only the coefficient changes -- so
    this reuses `core.project.frames` + `_torsion_grad` unmodified and stays O(n).
    """
    n = len(phi0)
    inv = 1.0 / np.asarray(sd, float)
    keep = np.ones(len(dhat), bool) if mask is None else np.asarray(mask, bool)

    def fg(x):
        phi = x[:n]; psi = x[n:]
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rvec = CA[i] - CA[j]
        d = np.maximum(np.sqrt((rvec * rvec).sum(1)), 1e-9)
        r = (d - dhat) * inv
        v, g = rho(r, kind, s)
        f = float(np.where(keep, v, 0.0).sum())
        #: d f / d d_ij = rho'(r) * dr/dd = rho'(r) / sd
        coef = (np.where(keep, g * inv, 0.0) / d)[:, None] * rvec
        gCA = np.zeros_like(CA)
        np.add.at(gCA, i, coef)
        np.add.at(gCA, j, -coef)
        return f, pj._torsion_grad(G, CA, gCA)

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    res = minimize(fg, x0, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return res.x[:n], res.x[n:], float(res.fun)


def fit_gnc(dhat, sd, i, j, phi0, psi0, kind="cauchy", schedule=(8.0, 4.0, 2.0, 1.0, 0.5),
            maxiter=200):
    """GRADUATED NON-CONVEXITY: anneal the loss scale from nearly-quadratic down to sharp.

    Why this arm exists. The 8-target smoke test showed the redescending `welsch` loss with the
    classic signature of its failure mode: a BETTER median (2.466 vs 2.777) and far more targets
    under 2 A (0.38 vs 0.12), with a WORSE mean (2.967 vs 2.704). A redescending loss solves many
    targets well and fails catastrophically on a few, and the mean is set by the failures.

    Two mechanisms produce that, and both are cured by annealing rather than by abandoning the
    loss. First, a redescending objective is non-convex, so a fit started far from the answer
    settles into a basin where it has already decided which restraints to ignore, and it never
    reconsiders. Second -- and this one is specific to our protocol -- the loss is BOUNDED ABOVE,
    so a structure that violates every restraint incurs only a bounded penalty, which corrupts
    the objective-only multi-start selection rule: a uniformly terrible fit can present a
    competitive objective value.

    Graduated non-convexity is the standard cure. Start at a scale so large the loss is
    effectively quadratic (every restraint pulls, the problem is near-convex, the basin is
    global), then anneal, warm-starting each stage from the previous solution, so restraints are
    discarded only once the structure is already approximately right. The final objective is
    evaluated at the sharpest scale, which restores a meaningful selection rule.
    """
    phi, psi = np.asarray(phi0, float), np.asarray(psi0, float)
    f = np.inf
    for s in schedule:
        phi, psi, f = fit_robust(dhat, sd, i, j, phi, psi, kind, s, None, maxiter)
    return phi, psi, f


def fit_trimmed(dhat, sd, i, j, phi0, psi0, q=0.10, maxiter=400):
    """Squared error, then DISCARD the worst q of restraints and refit. The NMR analogue.

    The trim is chosen by restraint violation, which is native-free -- no RMSD is consulted.
    """
    phi, psi, _ = fit_robust(dhat, sd, i, j, phi0, psi0, "squared", 1.0, None, maxiter)
    CA = np.asarray(I.build_ca(phi, psi), float)
    d = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    z = np.abs(d - dhat) / np.asarray(sd, float)
    cut = np.quantile(z, 1.0 - q)
    mask = z <= cut
    return fit_robust(dhat, sd, i, j, phi, psi, "squared", 1.0, mask, maxiter)


# --------------------------------------------------------------- leave-fold-out scale choice
def _residual_of_fit(t, dhat, sd, i, j, starts, kind, s, q=None):
    """The OBJECTIVE-selected fit's standardised restraint residual. Native-free."""
    best = None
    for phi0, psi0, _tag in starts:
        if kind == "trimmed":
            phi, psi, f = fit_trimmed(dhat, sd, i, j, phi0, psi0, q)
        else:
            phi, psi, f = fit_robust(dhat, sd, i, j, phi0, psi0, kind, s)
        if best is None or f < best[0]:
            best = (f, phi, psi)
    CA = np.asarray(I.build_ca(best[1], best[2]), float)
    d = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    #: score every arm on the SAME yardstick -- median |z| -- so scales are comparable
    return float(np.median(np.abs(d - dhat) / np.asarray(sd, float))), best[1], best[2]


def run(targets=None, n_start=6, use_debias=True):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    data = C.gather(tg)
    deb = {}
    if use_debias:
        for f in sorted(set(folds)):
            train = [p for p in pdbs if data[p]["fold"] != f]
            fn, _ = C.fit_correction(data, train, "sep")
            deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    # ---- per-target inputs, computed once
    prep = {}
    for t in tg:
        p = t["pdb"]; d = data[p]
        dhat = d["dhat"] - (deb[d["fold"]](d["sep"]) if use_debias else 0.0)
        rng = SD.stable_rng(p)
        prep[p] = {"dhat": np.maximum(dhat, 2.0), "sd": d["sd"], "i": d["i"], "j": d["j"],
                   "nat": d["nat"], "fold": d["fold"],
                   "starts": D.starts(p, d["seq"], d["n"], d["fold"], n_start, rng)}

    arms = ["squared"]
    for k in LOSSES[1:]:
        arms.append(k)
    arms += ["trimmed", "gnc_cauchy", "gnc_welsch", "ORACLE_scale"]

    res = {a: [] for a in arms}
    chosen = {}
    path = os.path.join(RESULTS, "robust.json")

    # ---- leave-fold-out scale selection, on the NATIVE-FREE residual only
    #: COST NOTE. The selection loop is 5 folds x 19 configurations x |train| targets x n_start
    #: fits. At |train| = 24 and n_start = 6 that is ~13,700 fits BEFORE any arm is evaluated,
    #: which on a loaded machine runs for hours. The selection only has to choose one scalar per
    #: loss per fold, and a scale chosen on 12 targets is not meaningfully worse than one chosen
    #: on 24 -- so the training slice and its start count are cut, and the ARMS are left at full
    #: strength. The reported arms still use `n_start` starts on all 126 targets.
    N_SEL, K_SEL = 12, 3
    print(f"selecting loss scales leave-fold-out on {N_SEL} training targets x {K_SEL} starts "
          f"(objective residual, never RMSD)", flush=True)
    sel = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if prep[p]["fold"] != f][:N_SEL]
        for k in LOSSES[1:]:
            best = None
            for s in SCALES:
                v = np.mean([_residual_of_fit(None, prep[p]["dhat"], prep[p]["sd"],
                                              prep[p]["i"], prep[p]["j"],
                                              prep[p]["starts"][:K_SEL], k, s)[0]
                             for p in train])
                if best is None or v < best[0]:
                    best = (v, s)
            sel[(f, k)] = best[1]
        best = None
        for q in TRIMS:
            v = np.mean([_residual_of_fit(None, prep[p]["dhat"], prep[p]["sd"],
                                          prep[p]["i"], prep[p]["j"],
                                          prep[p]["starts"][:K_SEL], "trimmed", 1.0, q)[0]
                         for p in train])
            if best is None or v < best[0]:
                best = (v, q)
        sel[(f, "trimmed")] = best[1]
        print(f"  fold {f}: " + ", ".join(f"{k}={sel[(f,k)]}" for k in
                                          list(LOSSES[1:]) + ["trimmed"]), flush=True)
    chosen = {f"{f}|{k}": sel[(f, k)] for (f, k) in sel}

    # ---- the arms
    for c, p in enumerate(pdbs):
        q = prep[p]
        f = q["fold"]
        for a in arms:
            if a == "ORACLE_scale":
                #: the in-fold ceiling: the best scale for THIS target, chosen on its own
                #: residual. ORACLE only in the sense of in-fold tuning; still no RMSD.
                cand = []
                for k in LOSSES[1:]:
                    for s in SCALES:
                        cand.append((_residual_of_fit(None, q["dhat"], q["sd"], q["i"],
                                                      q["j"], q["starts"][:K_SEL], k, s),
                                     k, s))
                v, k, s = min(cand, key=lambda z: z[0][0])
                phi, psi = v[1], v[2]
            elif a.startswith("gnc_"):
                best = None
                for phi0, psi0, _t in q["starts"]:
                    ph, ps, fv = fit_gnc(q["dhat"], q["sd"], q["i"], q["j"], phi0, psi0,
                                         kind=a[4:])
                    if best is None or fv < best[0]:
                        best = (fv, ph, ps)
                phi, psi = best[1], best[2]
            elif a == "trimmed":
                best = None
                for phi0, psi0, _t in q["starts"]:
                    ph, ps, fv = fit_trimmed(q["dhat"], q["sd"], q["i"], q["j"],
                                             phi0, psi0, sel[(f, "trimmed")])
                    if best is None or fv < best[0]:
                        best = (fv, ph, ps)
                phi, psi = best[1], best[2]
            else:
                s = 1.0 if a == "squared" else sel[(f, a)]
                best = None
                for phi0, psi0, _t in q["starts"]:
                    ph, ps, fv = fit_robust(q["dhat"], q["sd"], q["i"], q["j"],
                                            phi0, psi0, a, s)
                    if best is None or fv < best[0]:
                        best = (fv, ph, ps)
                phi, psi = best[1], best[2]
            res[a].append(float(I.ca_rmsd(I.build_ca(phi, psi), q["nat"])))
        if (c + 1) % 10 == 0 or c + 1 == len(pdbs):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1, "scales": chosen}, fh)
            print(f"  {c+1}/{len(pdbs)} checkpointed", flush=True)

    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray(res["squared"], float)
    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "scales": chosen, "arms": {},
           "per_target": {a: dict(zip(pdbs, map(float, res[a]))) for a in arms}}
    for a in arms:
        v = np.asarray(res[a], float)
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_squared": I.paired(v, base, folds=folds, names=pdbs),
                          "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_robust", out, n_expected=len(pdbs))

    print(f"\nincumbent {ref.mean():.3f}   (scales chosen leave-fold-out on the objective)\n")
    print(f"{'loss':<14}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs squared':>24}{'vs incumbent':>24}")
    for a in arms:
        s = out["arms"][a]
        vs, vi = s["vs_squared"], s["vs_incumbent"]
        print(f"{a:<14}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"  {vs['mean_diff']:+.3f} [{vs['ci95'][0]:+.3f},{vs['ci95'][1]:+.3f}]"
              f"  {vi['mean_diff']:+.3f} [{vi['ci95'][0]:+.3f},{vi['ci95'][1]:+.3f}]")
    print("\n(if no robust loss beats `squared`, the outlier reading of K7 is WRONG and the"
          "\n restraint error is diffuse rather than concentrated -- a clean negative worth "
          "having.)")
    return out


if __name__ == "__main__":
    run()
