"""s18/exp_leverage.py -- H8: does GEOMETRIC LEVERAGE re-weighting repair the refinement?

THE MOTIVATION, AND ITS RETRACTION -- BOTH RECORDED, BECAUSE THE ARM OUTLIVED ITS REASON.

The arm was pre-registered on the coordinator's `shuf_paired` control, which appeared to show that
the distogram's errors are harmful because of WHICH PAIRS they land on (2.072 A against the real
distogram's 3.610).  ADVERSARIAL then audited `s18/objceil.py` and found that `shuf_paired` passes
a PERMUTED WEIGHT VECTOR into the fit (`sd[pi]`, not `sd`), so it is not the deployed functional.
Against `shuffled` (2.609, residual permuted, weights kept in place) it differs only in whether
`sd` moves, so the extra 0.537 A is attributable to permuting the WEIGHTS, not to error
assignment.  **The coordinator has withdrawn that claim and the honest split is unmeasured.**

The arm is run anyway, at the coordinator's explicit request, exactly as pre-registered, so the
record shows an arm that outlived its motivation.  It is a legitimate native-free hypothesis on
its own terms.  **2.072 A is NOT its ceiling and is quoted nowhere as one**; the reference arms
that still stand are `shuffled` 2.609, `shuf_strat` 2.560, `isotropic` 2.573.

NOT MINE TO RUN.  ADVERSARIAL owns `wperm_only` (real dhat, permuted weights) and `wflat` (real
dhat, uniform weights).  This arm's own control permutes the LEVERAGE factor while leaving
`w = 1/sd^2` in place -- deliberately a different question -- so there is no duplication.

THE ARM.  Re-weight the refinement objective by how much each pair's TARGET distance can move the
emitted structure.  Native-free: everything is evaluated at the Control-A start.

    E(theta) = sum_p w_p (d_p(theta) - dhat_p)^2,        w_p = 1/sd_p^2

At the optimum `grad E = 0`.  Differentiating the stationarity condition in `dhat_p`,

    H dtheta*/ddhat_p = 2 w_p grad d_p,     H = 2 G^T diag(w) G   (Gauss-Newton, PSD)

so, with `J` the SUPERPOSED Ca Jacobian (rigid modes projected out, so the leverage is a real
structural displacement and not a frame artefact),

    lev_p = || J H^+ (2 w_p grad d_p) || / sqrt(n)      [Angstrom of structure per Angstrom of dhat_p]

and the arm minimises the SAME objective with

    w'_p = w_p * (lev_p / median_p lev)^(-2 beta),    beta in {0, 0.25, 0.5, 0.75, 1}

PRE-REGISTERED (`s18/PREREG_exp.md` H8) and NOT extended.  beta = 0 is exactly the deployed
weighting and must reproduce `refine_full`; that is the arm's validity gate.

CONTROLS, both mandatory, both run.
    lev_perm    the leverage vector PERMUTED ACROSS PAIRS -- identical weight distribution, the
                leverage->pair assignment destroyed.  This is the exact analogue of the
                coordinator's own shuffle and is what decides whether leverage carries
                information rather than merely re-scaling the objective.
    lev_rand    a random re-weighting matched to the realised w' distribution.

TWO CAUTIONS, adopted verbatim from the coordinator.
  1. Sprint 17 refuted 58 functionals of the distogram INCLUDING uncertainty weighting -- for
     SELECTION.  This is REFINEMENT: a different use and a different functional.  It inherits
     neither that refutation nor the alpha ceiling.
  2. A null here is a clean result, not a surprise -- and now more so than when the arm was
     written, since its motivating ceiling has been retracted as confounded.

AND ONE OF MY OWN.  There is no native-free rule that selects beta.  If an interior beta is best
that is a DIAGNOSTIC, not a deployable result, and it is labelled as one.

    python -m s18.exp_leverage
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                # noqa: E402
from s15 import align_lib as AL                # noqa: E402
from s15 import seed as SD                     # noqa: E402
from s18 import exp_run as XR                  # noqa: E402
from s18 import math_anova as MA               # noqa: E402

RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "exp_leverage.json")
BETAS = (0.0, 0.25, 0.5, 0.75, 1.0)            # PRE-REGISTERED.  Do not extend.


def leverage(phi0, psi0, i, j, w):
    """`(npair,)` structural leverage of each pair's TARGET distance, at the given structure.

    Returns Angstrom of Ca displacement per Angstrom of `dhat_p`, with the six rigid-body
    directions projected out of the Jacobian so the number is a real structural change.
    """
    n = len(phi0)
    #: `sup_jacobian` returns (J_sup (3n,2n) with the six rigid directions projected out,
    #: J_raw (n,3,2n), CA).  The SUPERPOSED one is what prices a real structural change; the
    #: raw one is what the pair-distance gradient needs, since a pair distance is itself
    #: rigid-invariant and must not be double-projected.
    Jf, Jr, CA = AL.sup_jacobian(phi0, psi0)
    rv = CA[i] - CA[j]
    d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
    uhat = rv / d[:, None]
    #: grad_theta d_p = uhat . (J_i - J_j)
    G = np.einsum("pa,pab->pb", uhat, Jr[i] - Jr[j])       # (npair, 2n)
    H = 2.0 * (G.T * w) @ G
    Hp = np.linalg.pinv(H, rcond=1e-10)
    dth = (2.0 * w)[:, None] * (G @ Hp)                    # (npair, 2n) = dtheta*/ddhat_p
    return np.linalg.norm(dth @ Jf.T, axis=1) / np.sqrt(n)


def run(targets=None, out=OUT):
    tg = targets if targets is not None else I.targets()
    deb = MA.debias_map(tg)
    rows = []
    if os.path.exists(out):
        try:
            p = json.load(open(out))
            if not p.get("complete"):
                rows = p["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()

    for t in tg:
        pdb, n = t["pdb"], int(t["n"])
        if pdb in done:
            continue
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        _W, _PH, _PS, avg, _ca, phi0, psi0 = XR.start_structure(t, nat)
        d = MA.gather_one(t, deb[int(t["fold"])])
        i, j, dhat, sd = d["i"], d["j"], d["dhat"], d["sd"]
        w = 1.0 / np.asarray(sd, float) ** 2
        rng = SD.stable_rng(pdb, "s18lever")

        lev = leverage(phi0, psi0, i, j, w)
        rel = lev / max(np.median(lev), 1e-12)
        e = {"pdb": pdb, "n": n, "fold": int(t["fold"]),
             "avg": float(I.ca_rmsd(avg, nat)),
             "lev_median": float(np.median(lev)), "lev_iqr": float(np.subtract(
                 *np.percentile(lev, [75, 25]))), "lev_max": float(lev.max()),
             "lev_cv": float(lev.std() / max(lev.mean(), 1e-12))}
        #: ORACLE DIAGNOSTIC (never a result): is the distogram's error actually correlated
        #: with leverage?  This is the mechanism the coordinator's shuf_paired implies.
        dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
        res = np.abs(np.asarray(dhat, float) - dtrue)
        rk = lambda x: np.argsort(np.argsort(x)).astype(float)      # noqa: E731
        e["ORC_rho_lev_absresid"] = float(np.corrcoef(rk(lev), rk(res))[0, 1])
        e["ORC_rho_w_absresid"] = float(np.corrcoef(rk(w), rk(res))[0, 1])
        e["ORC_rho_lev_sep"] = float(np.corrcoef(rk(lev), rk((j - i).astype(float)))[0, 1])

        for b in BETAS:
            k = f"b{int(round(b * 100)):03d}"
            wp = rel ** (-2.0 * b)
            p_, q_, f_ = AL.fit(dhat, sd, i, j, phi0, psi0, wpair=wp)
            e[k] = float(I.ca_rmsd(I.build_ca(p_, q_), nat))
            e[k + "_obj"] = float(((((np.sqrt((
                (I.build_ca(p_, q_)[i] - I.build_ca(p_, q_)[j]) ** 2).sum(1))
                - dhat) / sd) ** 2)).sum())
            #: CONTROL 1 -- identical weight DISTRIBUTION, leverage->pair assignment destroyed
            wq = wp[rng.permutation(len(wp))]
            p_, q_, _f = AL.fit(dhat, sd, i, j, phi0, psi0, wpair=wq)
            e[k + "_perm"] = float(I.ca_rmsd(I.build_ca(p_, q_), nat))
            #: CONTROL 2 -- matched-random weights with the same mean and sd in log space
            lg = np.log(np.maximum(wp, 1e-12))
            wr = np.exp(rng.standard_normal(len(wp)) * lg.std() + lg.mean())
            p_, q_, _f = AL.fit(dhat, sd, i, j, phi0, psi0, wpair=wr)
            e[k + "_rand"] = float(I.ca_rmsd(I.build_ca(p_, q_), nat))

        rows.append(e)
        json.dump({"rows": rows, "complete": False, "betas": list(BETAS)}, open(out, "w"))
        if len(rows) % 10 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s) {pdb}", flush=True)

    json.dump({"rows": rows, "complete": len(rows) == len(tg), "betas": list(BETAS),
               "n_expected": len(tg), "n_rows": len(rows)}, open(out, "w"))
    print(f"DONE leverage {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


def report(path=OUT):
    from s18.exp_report import boot_fold
    o = json.load(open(path))
    if not o.get("complete"):
        print(f"** exp_leverage.json PARTIAL: {len(o['rows'])} rows -- NOT final **")
    rows = o["rows"]
    rng = SD.stable_rng("s18", "leverage_report")
    folds = np.array([r["fold"] for r in rows])
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    n = len(rows)
    print("=" * 100)
    print(f"[H8] GEOMETRIC-LEVERAGE RE-WEIGHTING, n = {n}.  beta PRE-REGISTERED, not extended.")
    print("     w'_p = w_p * (lev_p / median lev)^(-2 beta).  beta = 0 IS the deployed weighting.")
    print("=" * 100)
    print(f"\n  {'beta':>6}{'RMSD':>9}{'median':>9}{'vs avg [95% CI]':>28}{'W/L':>9}"
          f"{'permCTRL':>10}{'randCTRL':>10}{'vs perm':>9}")
    for b in o["betas"]:
        k = f"b{int(round(b * 100)):03d}"
        v = g(k)
        mu, lo, hi = boot_fold(v - g("avg"), folds, rng)
        dp = v - g(k + "_perm")
        print(f"  {b:>6.2f}{v.mean():>9.3f}{np.median(v):>9.3f}"
              f"   {mu:+.3f} [{lo:+.3f},{hi:+.3f}]"
              f"{int((v < g('avg')).sum()):>5}/{int((v > g('avg')).sum()):<4}"
              f"{g(k + '_perm').mean():>10.3f}{g(k + '_rand').mean():>10.3f}"
              f"{dp.mean():>+9.3f}")
    base = g("b000")
    print(f"\n  VALIDITY GATE: beta = 0 must equal `refine_full`.  mean {base.mean():.3f}")
    s17 = os.path.join(ROOT, "s17", "results", "refine.json")
    if os.path.exists(s17):
        m = {r["pdb"]: r for r in json.load(open(s17))["rows"]}
        com = [r for r in rows if r["pdb"] in m]
        if com:
            dd = np.array([r["b000"] - m[r["pdb"]]["refine_full"] for r in com])
            print(f"     vs s17 refine_full on {len(com)} shared: mean|d| {np.abs(dd).mean():.5f}"
                  f"  -> {'PASS' if np.abs(dd).mean() < 0.01 else 'FAIL'}")
    for b in o["betas"][1:]:
        k = f"b{int(round(b * 100)):03d}"
        mu, lo, hi = boot_fold(g(k) - base, folds, rng)
        mp, lp, hp = boot_fold(g(k) - g(k + "_perm"), folds, rng)
        print(f"  beta={b}: vs beta=0 {mu:+.3f} [{lo:+.3f},{hi:+.3f}]   "
              f"vs its PERMUTED control {mp:+.3f} [{lp:+.3f},{hp:+.3f}]")
    print("\n  ORACLE DIAGNOSTIC -- is the distogram's error actually where leverage is?")
    print(f"     rho(leverage, |residual|)      {g('ORC_rho_lev_absresid').mean():+.3f}"
          f"   median {np.median(g('ORC_rho_lev_absresid')):+.3f}")
    print(f"     rho(1/sd^2,   |residual|)      {g('ORC_rho_w_absresid').mean():+.3f}")
    print(f"     rho(leverage, separation)      {g('ORC_rho_lev_sep').mean():+.3f}")
    print(f"     leverage spread: CV {g('lev_cv').mean():.3f}, max/median "
          f"{(g('lev_max') / np.maximum(g('lev_median'), 1e-12)).mean():.1f}")
    print("\n  CAUTIONS carried from the coordinator, and they are not softened:")
    print("   * Sprint 17 refuted 58 distogram functionals INCLUDING uncertainty weighting --")
    print("     for SELECTION.  This is REFINEMENT.  It inherits neither that refutation nor the")
    print("     alpha ceiling, and is read on its own evidence.")
    print("   * The motivating ceiling (shuf_paired, 2.072 A) has been RETRACTED as confounded:")
    print("     it passed a permuted WEIGHT vector into the fit, so it is not the deployed")
    print("     functional.  This arm outlived its motivation and is reported on its own")
    print("     evidence.  A null here is a clean result.")
    print("   * No native-free rule selects beta.  An interior best is a DIAGNOSTIC, not a")
    print("     deployable result.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
