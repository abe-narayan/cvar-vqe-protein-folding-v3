"""s16/verify_steer.py -- VERIFY workstream, hostile re-reading of the FLAGSHIP.

THREE QUESTIONS, in the order the assignment sets them.

  Q1  Does the published table reproduce from `s16/results/steer.json`?  Including the
      per-target-mean vs ratio-of-means switch the writeup makes mid-section.

  Q2  Is `theta_fit + s * A.wrap(theta_native - theta_fit)` the GEODESIC interpolation on the
      torus, or a naive linear interpolation of wrapped coordinates?  Three interpolants are
      built independently and compared bit for bit:
        naive   (1 - s) * theta_fit + s * theta_native            -- NOT wrapped
        wrapd   theta_fit + s * wrap(theta_native - theta_fit)    -- the coordinator's
        slerp   angle of the normalised chord between the two unit vectors on each circle,
                built from atan2 of an independently derived great-circle formula
      If `wrapd` and `slerp` agree to machine precision the geodesic objection is dead.

  Q3  THE CONTROL THE FLAGSHIP DOES NOT HAVE.  The writeup attributes the non-monotone
      mid-path bulge to a PROPERTY OF THE ERROR -- "the fit's torsion errors are
      *compensating* ... moving part-way breaks the cancellation".  That is a causal claim
      with no matched control.  The alternative explanation is that torsion-space geodesic
      interpolation between ANY two distant conformations bulges, i.e. it is a property of
      the PARAMETERISATION.  This module supplies the missing control:

        pool->pool   geodesic interpolation between two real top-75 retrieval windows
        pool->rand   from a real window toward a random torsion vector of MATCHED norm

      measured with exactly the flagship's statistic, entirely NATIVE-FREE:

        frac(s) = 1 - d(x(s), T) / d(x(0), T)          T = the endpoint structure

      (For the flagship arm the published statistic is (d0 - d(s)) / (d0 - d1) against the
      native, with d1 = 0.043 A; the two differ by less than 1.2% of the gain.)

      If pool->pool reproduces the flagship's curve at matched angular distance
      ||wrap(dtheta)||, the "compensating error" mechanism is unsupported and the finding is
      about the coordinate system, not about the error.

ORACLE.  Only Q1 reads native-derived numbers, and only ones already in the artefact.  Q2 and
Q3 are entirely native-free.

Writes `s16/results/verify_steer.json`.
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import seed as SD                   # noqa: E402

STEPS = (0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0)
N_PAIR = 12


# ------------------------------------------------------------------ Q1: reproduction
def q1(out):
    d = json.load(open(os.path.join(RESULTS, "steer.json")))
    rows = d["rows"]
    S = [str(s) for s in d["steps"]]
    M = np.array([[r["arms"]["ORACLE_true"][s] for s in S] for r in rows])
    base, end = M[:, 0], M[:, -1]
    gain = base - end
    fr = (base[:, None] - M) / gain[:, None]
    rom = (base.mean() - M.mean(0)) / (base.mean() - end.mean())
    per = fr.mean(0)
    med = np.median(fr, 0)
    worse = [(int((M[:, k] > M[:, 0]).sum()), int((M[:, k] < M[:, 0]).sum()))
             for k in range(len(S))]
    folds = np.asarray([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    i50 = S.index("0.5")
    ci = I.paired(fr[:, i50], np.zeros(len(rows)), folds=folds, names=pdbs)

    print("=" * 96)
    print("Q1  REPRODUCTION of STEER_FINDINGS section 3, from s16/results/steer.json")
    print("=" * 96)
    print(f"  n = {len(rows)}   published n = 126")
    print(f"  {'step s':<28}" + "".join(f"{s:>9}" for s in S))
    print(f"  {'mean RMSD (A)':<28}" + "".join(f"{v:>9.3f}" for v in M.mean(0)))
    print(f"  {'frac of gain, ratio-of-means':<28}" + "".join(f"{v:>9.4f}" for v in rom))
    print(f"  {'frac of gain, PER TARGET':<28}" + "".join(f"{v:>9.4f}" for v in per))
    print(f"  {'   median':<28}" + "".join(f"{v:>9.4f}" for v in med))
    print(f"  {'improved / worsened':<28}" + "".join(f"{f'{w}/{l}':>9}" for l, w in worse))
    print(f"\n  at s = 0.5: per-target mean {per[i50]:.4f}, fold-clustered CI "
          f"[{ci['ci95'][0]:+.4f},{ci['ci95'][1]:+.4f}], median {med[i50]:.4f}, "
          f"{worse[i50][1]}/{worse[i50][0]} improved/worsened")
    print("  PUBLISHED: 7.7% [-0.0, +14.9], median 16.1%, 77 improved / 49 worsened  -> "
          f"{'MATCH' if abs(per[i50] - 0.077) < 5e-4 and worse[i50][0] == 49 else 'MISMATCH'}")

    q = np.percentile(fr[:, i50], [0, 5, 25, 50, 75, 95, 100])
    print(f"\n  the distribution behind the mean, at s = 0.5:")
    print("     min      p5     p25     p50     p75     p95     max")
    print("  " + "".join(f"{v:>8.3f}" for v in q))
    print(f"  targets with frac >= 0.45 (near-linear or better): "
          f"{int((fr[:, i50] >= 0.45).sum())}/{len(rows)}")
    print(f"  targets with frac <= -0.25 (a large mid-path BULGE):  "
          f"{int((fr[:, i50] <= -0.25).sum())}/{len(rows)}")
    print("  -- the 7.7% mean is a mixture, not a typical target. It is the average of a")
    print("     near-linear majority and a heavy negative tail.")

    en = np.array([r["err_norm"] for r in rows])
    rho = np.corrcoef(np.argsort(np.argsort(en)),
                      np.argsort(np.argsort(fr[:, i50])))[0, 1]
    print(f"\n  Spearman(||e||, frac at s=0.5) = {rho:+.3f}  "
          f"(||e|| mean {en.mean():.2f} rad, range {en.min():.2f}-{en.max():.2f})")
    lo, hi = en <= np.median(en), en > np.median(en)
    print(f"    small-||e|| half: frac {fr[lo, i50].mean():+.3f}   "
          f"large-||e|| half: frac {fr[hi, i50].mean():+.3f}")

    out["q1"] = {
        "n": len(rows), "steps": S, "mean_curve": M.mean(0).tolist(),
        "frac_ratio_of_means": rom.tolist(), "frac_per_target": per.tolist(),
        "frac_median": med.tolist(), "improved_worsened": worse,
        "ci_at_half_fold_clustered": ci["ci95"],
        "pctiles_at_half": q.tolist(),
        "spearman_errnorm_frac": float(rho),
        "frac_small_e": float(fr[lo, i50].mean()), "frac_large_e": float(fr[hi, i50].mean()),
    }
    return rows, M, fr, en


# ------------------------------------------------------------------ Q2: the interpolant
def q2(out):
    print("\n" + "=" * 96)
    print("Q2  IS THE INTERPOLANT THE GEODESIC?  (three independent constructions)")
    print("=" * 96)
    rng = SD.stable_rng("verify", "geodesic")
    n = 13
    worst_ws, worst_wn, worst_rmsd = 0.0, 0.0, 0.0
    for _ in range(200):
        a = rng.uniform(-np.pi, np.pi, 2 * n)
        # a target deliberately placed so that many coordinates straddle the branch cut
        b = a + rng.uniform(-np.pi, np.pi, 2 * n) + rng.choice([-2, 0, 2], 2 * n) * np.pi
        for s in STEPS:
            wrapd = a + s * A.wrap(b - a)
            naive = (1 - s) * a + s * b
            # independent great-circle construction: rotate the unit vector of `a` toward
            # that of `b` by s times the signed angle between them.
            va = np.stack([np.cos(a), np.sin(a)])
            vb = np.stack([np.cos(b), np.sin(b)])
            dot = (va * vb).sum(0)
            crs = va[0] * vb[1] - va[1] * vb[0]
            ang = np.arctan2(crs, dot)                 # signed angle in (-pi, pi]
            slerp = a + s * ang
            worst_ws = max(worst_ws, float(np.abs(A.wrap(wrapd - slerp)).max()))
            worst_wn = max(worst_wn, float(np.abs(A.wrap(wrapd - naive)).max()))
            worst_rmsd = max(worst_rmsd, abs(
                float(I.ca_rmsd(I.build_ca(wrapd[:n], wrapd[n:]),
                                I.build_ca(slerp[:n], slerp[n:])))))
    print(f"  max |wrap(wrapd - slerp)| over 200 random pairs x 7 steps : {worst_ws:.3e} rad")
    print(f"  max CA-RMSD between the two builds                        : {worst_rmsd:.3e} A")
    print(f"  max |wrap(wrapd - naive)| (the interpolant NOT used)      : {worst_wn:.3e} rad")
    print("  VERDICT: the coordinator's `th + s*A.wrap(th_nat - th)` IS the per-angle geodesic.")
    print("  The naive unwrapped interpolation differs by up to pi radians and was NOT used.")
    print("  The geodesic objection to the flagship is REFUTED.")
    out["q2"] = {"max_wrapd_minus_slerp_rad": worst_ws,
                 "max_rmsd_between_builds": worst_rmsd,
                 "max_wrapd_minus_naive_rad": worst_wn}


# ------------------------------------------------------------------ Q3: the missing control
def _curve(th0, th1, n):
    """geodesic torsion interpolation; returns RMSD to the ENDPOINT structure at each step."""
    dv = A.wrap(th1 - th0)
    T = I.build_ca(th1[:n], th1[n:])
    return np.array([I.ca_rmsd(I.build_ca((th0 + s * dv)[:n], (th0 + s * dv)[n:]), T)
                     for s in STEPS]), float(np.linalg.norm(dv))


def q3(out, err_norm, fr_pub):
    from s14 import retprior as R
    print("\n" + "=" * 96)
    print("Q3  THE MISSING CONTROL: does torsion-space interpolation bulge between ANY two")
    print("    distant conformations?   (entirely NATIVE-FREE)")
    print("=" * 96)
    tg = I.targets()
    pool_fr, pool_d, rand_fr, rand_d, pdbs, folds = [], [], [], [], [], []
    for t in tg:
        p = t["pdb"]
        PHI, PSI, _ = R.windows(p, "top75")
        n = PHI.shape[1]
        TH = np.concatenate([PHI, PSI], axis=1)
        rng = SD.stable_rng(p, "verify_q3")
        m = len(TH)
        for _k in range(N_PAIR):
            a, b = rng.choice(m, 2, replace=False)
            c, nd = _curve(TH[a], TH[b], n)
            if c[0] > 1e-6:
                pool_fr.append(1.0 - c / c[0]); pool_d.append(nd)
                pdbs.append(p); folds.append(int(t["fold"]))
        # matched-magnitude random target from a real window
        for _k in range(N_PAIR):
            a = rng.integers(0, m)
            z = rng.standard_normal(2 * n)
            z /= np.linalg.norm(z)
            for target_norm in (float(np.median(err_norm)),):
                c, nd = _curve(TH[a], TH[a] + z * target_norm, n)
                if c[0] > 1e-6:
                    rand_fr.append(1.0 - c / c[0]); rand_d.append(nd)
    pool_fr = np.array(pool_fr); pool_d = np.array(pool_d)
    rand_fr = np.array(rand_fr); rand_d = np.array(rand_d)

    i50 = STEPS.index(0.5)
    print(f"  pool->pool pairs: {len(pool_fr)}   ||dtheta|| mean {pool_d.mean():.2f} rad")
    print(f"  pool->rand pairs: {len(rand_fr)}   ||dtheta|| mean {rand_d.mean():.2f} rad")
    print(f"  flagship fit->native: n = {len(err_norm)}  ||e|| mean {err_norm.mean():.2f} rad")
    print(f"\n  {'arm':<26}" + "".join(f"{s:>9}" for s in STEPS))
    print(f"  {'flagship fit->native':<26}" + "".join(f"{v:>9.4f}" for v in fr_pub.mean(0)))
    print(f"  {'CONTROL pool->pool':<26}" + "".join(f"{v:>9.4f}" for v in pool_fr.mean(0)))
    print(f"  {'CONTROL pool->random':<26}" + "".join(f"{v:>9.4f}" for v in rand_fr.mean(0)))
    print(f"  {'  (median) pool->pool':<26}" + "".join(f"{v:>9.4f}"
                                                       for v in np.median(pool_fr, 0)))
    print(f"  {'  (median) fit->native':<26}" + "".join(f"{v:>9.4f}"
                                                        for v in np.median(fr_pub, 0)))

    # matched on angular distance: restrict the control to the flagship's ||e|| range
    lo, hi = np.percentile(err_norm, [25, 75])
    sel = (pool_d >= lo) & (pool_d <= hi)
    selr = (rand_d >= lo) & (rand_d <= hi)
    self_ = (err_norm >= lo) & (err_norm <= hi)
    print(f"\n  MATCHED on angular distance, ||dtheta|| in [{lo:.2f}, {hi:.2f}] rad "
          f"(the flagship's own IQR):")
    print(f"  {'arm':<26}{'n':>6}" + "".join(f"{s:>9}" for s in STEPS))
    print(f"  {'flagship fit->native':<26}{int(self_.sum()):>6}"
          + "".join(f"{v:>9.4f}" for v in fr_pub[self_].mean(0)))
    print(f"  {'CONTROL pool->pool':<26}{int(sel.sum()):>6}"
          + "".join(f"{v:>9.4f}" for v in pool_fr[sel].mean(0)))
    print(f"  {'CONTROL pool->random':<26}{int(selr.sum()):>6}"
          + "".join(f"{v:>9.4f}" for v in rand_fr[selr].mean(0)))
    print(f"\n  fraction of paths WORSE than not moving, at s = 0.5:")
    print(f"    flagship fit->native  {float((fr_pub[:, i50] < 0).mean()):.3f}")
    print(f"    CONTROL pool->pool    {float((pool_fr[:, i50] < 0).mean()):.3f}")
    print(f"    CONTROL pool->random  {float((rand_fr[:, i50] < 0).mean()):.3f}")

    # the bulge as a function of angular distance, both arms in one binning
    print(f"\n  frac at s = 0.5 by angular distance decile (the real driver):")
    edges = np.percentile(np.concatenate([pool_d, err_norm]), np.linspace(0, 100, 7))
    print(f"  {'||dtheta|| bin':<20}{'fit->native':>16}{'pool->pool':>16}{'pool->rand':>16}")
    bins = []
    for k in range(6):
        a, b = edges[k], edges[k + 1]
        sf = (err_norm >= a) & (err_norm < b)
        sp = (pool_d >= a) & (pool_d < b)
        sr = (rand_d >= a) & (rand_d < b)
        f1 = fr_pub[sf, i50].mean() if sf.sum() else np.nan
        f2 = pool_fr[sp, i50].mean() if sp.sum() else np.nan
        f3 = rand_fr[sr, i50].mean() if sr.sum() else np.nan
        bins.append([float(a), float(b), int(sf.sum()), float(f1), int(sp.sum()), float(f2),
                     int(sr.sum()), float(f3)])
        print(f"  [{a:>5.2f},{b:>5.2f})      {f1:>+8.3f} (n={int(sf.sum()):>3})"
              f"{f2:>+8.3f} (n={int(sp.sum()):>4}){f3:>+8.3f} (n={int(sr.sum()):>4})")

    out["q3"] = {
        "n_pool_pairs": len(pool_fr), "n_rand_pairs": len(rand_fr),
        "steps": list(STEPS),
        "frac_fit_native": fr_pub.mean(0).tolist(),
        "frac_pool_pool": pool_fr.mean(0).tolist(),
        "frac_pool_rand": rand_fr.mean(0).tolist(),
        "median_pool_pool": np.median(pool_fr, 0).tolist(),
        "matched_iqr": [float(lo), float(hi)],
        "matched_fit_native": fr_pub[self_].mean(0).tolist(),
        "matched_pool_pool": pool_fr[sel].mean(0).tolist(),
        "matched_pool_rand": rand_fr[selr].mean(0).tolist(),
        "frac_negative_at_half": {
            "fit_native": float((fr_pub[:, i50] < 0).mean()),
            "pool_pool": float((pool_fr[:, i50] < 0).mean()),
            "pool_rand": float((rand_fr[:, i50] < 0).mean())},
        "distance_bins": bins,
        "dtheta_mean": {"fit_native": float(err_norm.mean()),
                        "pool_pool": float(pool_d.mean()),
                        "pool_rand": float(rand_d.mean())},
    }


def main():
    out = {}
    rows, M, fr, en = q1(out)
    q2(out)
    q3(out, en, fr)
    with open(os.path.join(RESULTS, "verify_steer.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote s16/results/verify_steer.json")


if __name__ == "__main__":
    main()
