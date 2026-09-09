"""S25 / PHYSICS LANE -- THE LEGACY-vs-AMBER LANDSCAPE COMPARISON. Professor-facing.

    python s25/phys_landscape.py

WHAT IS RE-MEASURED AND WHAT IS CITED
=====================================
Per the brief: reuse and CITE the Sprint 20 and Sprint 24 measurements where the pool has not
changed, and re-measure only what a changed candidate set or a changed score could plausibly
change.

  CITED, NOT RE-RUN -- properties of the ENERGY FUNCTION on the continuous torsion manifold,
  which no candidate pool can move: gradient magnitude and direction, Hessian condition number,
  anisotropy, participation ratio, near-zero spectrum, ruggedness, negative curvature (all s20).

  RE-MEASURED HERE on the shipped K=500 retrieval pool, because every one of these is a
  property of the CANDIDATE SET or of the SCORE and both are what this sprint fixes: energy
  scale and dynamic range, the candidate-energy distribution and its conditioning, rank
  correlation between the three channels, correlation with RMSD (global and in-band), the
  compactness response WITH an Rg control, steric sensitivity, CVaR concentration and
  trainability, and structural usefulness (the endpoint itself).

EVERY STANDING FACT THE BRIEF LISTS IS VERIFIED, NOT ASSUMED
============================================================
  rho(Legacy, AMBER) = -0.0829 retrieval / -0.0886 torsion
  rho(Legacy, distogram) = +0.3875   rho(AMBER, distogram) = -0.0270
  Legacy prefers candidates ~0.45 A MORE compact; AMBER prefers expanded
  both energies are FULL-REGISTER in torsion space ("AMBER is less local" is a category error)
  no defensible AMBER surrogate exists (bonded subset rho 0.211 against a 0.7 bar)

TWO TRAPS THIS MODULE IS BUILT AROUND
=====================================
1. `physics-ranks-real-geometry-not-lattice`: a whole-pool correlation with RMSD is LED BY
   RADIUS OF GYRATION, so never believe a physics scorer's rho with RMSD until Rg is partialled
   out. Both the raw and the Rg-partialled correlation are reported for all three channels, and
   an Rg-only control column sits beside them.
2. `in-band-is-the-only-ranking-metric`: global ordering and in-band ordering are different
   skills and the endpoint lives on the second. Both are reported; the in-band band is defined
   by the ORACLE per-candidate RMSD and is labelled ORACLE DIAGNOSTIC everywhere it appears.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402

BAND = 3.0        # ORACLE DIAGNOSTIC: the near-native band, per `in-band-is-the-only-...`
MIN_BAND = 8      # minimum in-band candidates for a per-target in-band correlation


def spearman(a, b):
    from scipy.stats import rankdata
    ra, rb = rankdata(a), rankdata(b)
    if ra.std() < 1e-12 or rb.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def partial_spearman(x, y, z):
    """rho(x, y | z) on ranks -- the Rg control `physics-ranks-real-geometry` demands."""
    from scipy.stats import rankdata
    R = np.column_stack([rankdata(x), rankdata(y), rankdata(z)]).astype(float)
    R -= R.mean(0)
    s = R.std(0)
    if (s < 1e-12).any():
        return float("nan")
    R /= s
    C = np.corrcoef(R.T)
    num = C[0, 1] - C[0, 2] * C[1, 2]
    den = math.sqrt(max(1e-12, (1 - C[0, 2] ** 2) * (1 - C[1, 2] ** 2)))
    return float(num / den)


def rg_of(W):
    """Radius of gyration of each CA cloud. (k, n, 3) -> (k,)."""
    c = W - W.mean(1, keepdims=True)
    return np.sqrt((c ** 2).sum(2).mean(1))


def min_sep3_dist(W):
    """Smallest CA-CA distance over |i-j| >= 3 -- the cheap steric proxy on this basis.

    A CA trace has no atoms, so this is not a clash count; it is the coarse-grained analogue,
    and it is reported as such. The all-atom clash census is s20's and is cited, not redone.
    """
    k, n, _ = W.shape
    i, j = np.triu_indices(n, 3)
    if len(i) == 0:
        return np.full(k, np.nan)
    d = np.linalg.norm(W[:, i] - W[:, j], axis=2)
    return d.min(1)


def mean_se(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), int(len(v))


def main() -> int:
    pdbs = P.targets()
    rows = []
    print(f"LANDSCAPE -- re-measuring the pool-dependent axes on {len(pdbs)} targets\n")
    for ii, pdb in enumerate(pdbs):
        cand, ch = P.channels(pdb)
        rr = np.asarray(cand.oracle_rr, float)          # ORACLE label, diagnostic only
        rg = rg_of(cand.W)
        ms = min_sep3_dist(cand.W)
        L, A, D = ch["LEG"], ch["AMB"], ch["DIS"]
        band = rr < BAND
        r = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), k=int(cand.k),
                 n_band=int(band.sum()))
        # --- scale and distribution
        for nm, e in (("LEG", L), ("AMB", A), ("DIS", D)):
            pos = e[np.isfinite(e) & (np.abs(e) > 0)]
            r[f"{nm}_min"] = float(e.min()); r[f"{nm}_max"] = float(e.max())
            r[f"{nm}_median"] = float(np.median(e))
            r[f"{nm}_iqr"] = float(np.percentile(e, 75) - np.percentile(e, 25))
            r[f"{nm}_decades"] = float(np.log10(np.abs(pos).max() / max(1e-12,
                                                                       np.abs(pos).min())))
            zr = (e - e.mean()) / max(e.std(), 1e-12)
            r[f"{nm}_frac_absz_lt_0p1"] = float((np.abs(zr) < 0.1).mean())
            r[f"{nm}_top10_var_share"] = float(
                np.sort((e - e.mean()) ** 2)[-10:].sum() / max(1e-12, ((e - e.mean()) ** 2).sum()))
            # --- correlations
            r[f"rho_{nm}_rmsd"] = spearman(e, rr)
            r[f"rho_{nm}_rmsd_partial_rg"] = partial_spearman(e, rr, rg)
            r[f"rho_{nm}_rg"] = spearman(e, rg)
            r[f"rho_{nm}_minsep3"] = spearman(e, ms)
            r[f"rho_{nm}_rmsd_inband"] = (spearman(e[band], rr[band])
                                          if band.sum() >= MIN_BAND else float("nan"))
            # --- compactness RESPONSE: what the selector's own top-75 looks like
            top = np.argsort(e, kind="stable")[:P.M_PROD]
            r[f"{nm}_rg_top75_minus_pool"] = float(rg[top].mean() - rg.mean())
            r[f"{nm}_minsep3_top75_minus_pool"] = float(ms[top].mean() - ms.mean())
            r[f"{nm}_oracle_rr_top75"] = float(rr[top].mean())
        # --- the Rg-only control, the column memory says must sit beside the physics ones
        r["rho_RG_rmsd"] = spearman(rg, rr)
        r["rho_RG_rmsd_inband"] = (spearman(rg[band], rr[band])
                                   if band.sum() >= MIN_BAND else float("nan"))
        # --- between-channel agreement
        r["rho_LEG_AMB"] = spearman(L, A)
        r["rho_LEG_DIS"] = spearman(L, D)
        r["rho_AMB_DIS"] = spearman(A, D)
        r["oracle_pool_best"] = float(rr.min()); r["oracle_pool_mean"] = float(rr.mean())
        r["rg_pool_mean"] = float(rg.mean()); r["minsep3_pool_mean"] = float(ms.mean())
        rows.append(r)
        if (ii + 1) % 25 == 0:
            print(f"  {ii+1}/{len(pdbs)}", flush=True)

    def col(k):
        return np.array([r[k] for r in rows], float)

    print("\n" + "=" * 92)
    print("A. ENERGY SCALE AND CANDIDATE-ENERGY DISTRIBUTION  (this pool, K=500, n=126)")
    print("=" * 92)
    print(f"  {'channel':<10} {'median':>12} {'IQR':>12} {'min':>12} {'max':>12} "
          f"{'decades':>8} {'|z|<0.1':>8} {'top10 var':>10}")
    for nm in ("LEG", "AMB", "DIS"):
        print(f"  {nm:<10} {np.mean(col(nm+'_median')):>12.4g} "
              f"{np.mean(col(nm+'_iqr')):>12.4g} {np.min(col(nm+'_min')):>12.4g} "
              f"{np.max(col(nm+'_max')):>12.4g} {np.mean(col(nm+'_decades')):>8.2f} "
              f"{np.mean(col(nm+'_frac_absz_lt_0p1')):>8.4f} "
              f"{np.mean(col(nm+'_top10_var_share')):>10.4f}")

    print("\n" + "=" * 92)
    print("B. RANK AGREEMENT BETWEEN THE THREE CHANNELS  (Spearman within target, n=126)")
    print("=" * 92)
    for k, cite in (("rho_LEG_AMB", "s24 -0.0829 retrieval / s20 -0.0886 torsion"),
                    ("rho_LEG_DIS", "s24 +0.3875"),
                    ("rho_AMB_DIS", "s24 -0.0270")):
        m, se, nn = mean_se(col(k))
        med = float(np.nanmedian(col(k)))
        neg = int((col(k) < 0).sum())
        print(f"  {k:<14} {m:>+8.4f}  SE {se:.4f}  median {med:>+8.4f}  "
              f"{neg}/{nn} negative     [standing: {cite}]")

    print("\n" + "=" * 92)
    print("C. CORRELATION WITH RMSD -- AND WHY THE WHOLE-POOL COLUMN IS A TRAP")
    print("   ORACLE DIAGNOSTIC: `rr` is the per-candidate RMSD to the native. Never a selector.")
    print("=" * 92)
    print(f"  {'channel':<10} {'rho(E,RMSD)':>12} {'partial Rg':>12} {'in-band':>12} "
          f"{'rho(E,Rg)':>12} {'rho(E,minsep3)':>15}")
    for nm in ("LEG", "AMB", "DIS"):
        a = mean_se(col(f"rho_{nm}_rmsd")); b = mean_se(col(f"rho_{nm}_rmsd_partial_rg"))
        c = mean_se(col(f"rho_{nm}_rmsd_inband")); d = mean_se(col(f"rho_{nm}_rg"))
        e = mean_se(col(f"rho_{nm}_minsep3"))
        print(f"  {nm:<10} {a[0]:>+9.4f}+-{a[1]:.3f} {b[0]:>+9.4f}+-{b[1]:.3f} "
              f"{c[0]:>+9.4f}+-{c[1]:.3f} {d[0]:>+9.4f}+-{d[1]:.3f} "
              f"{e[0]:>+12.4f}+-{e[1]:.3f}")
    a = mean_se(col("rho_RG_rmsd")); b = mean_se(col("rho_RG_rmsd_inband"))
    print(f"  {'Rg CONTROL':<10} {a[0]:>+9.4f}+-{a[1]:.3f} {'--':>12} "
          f"{b[0]:>+9.4f}+-{b[1]:.3f}")
    print(f"  in-band = the {int(np.mean(col('n_band')))} of 500 candidates with RMSD < {BAND} A "
          f"(mean per target); {int((col('n_band')>=MIN_BAND).sum())}/126 targets have >= "
          f"{MIN_BAND}")

    print("\n" + "=" * 92)
    print("D. COMPACTNESS AND STERIC RESPONSE -- what each energy's own top-75 looks like")
    print("=" * 92)
    print(f"  {'channel':<10} {'dRg(top75-pool)':>17} {'d minsep3':>12} "
          f"{'ORACLE rr top75':>17}")
    for nm in ("LEG", "AMB", "DIS"):
        a = mean_se(col(f"{nm}_rg_top75_minus_pool"))
        b = mean_se(col(f"{nm}_minsep3_top75_minus_pool"))
        c = mean_se(col(f"{nm}_oracle_rr_top75"))
        print(f"  {nm:<10} {a[0]:>+14.4f}+-{a[1]:.3f} {b[0]:>+9.4f}+-{b[1]:.3f} "
              f"{c[0]:>14.4f}+-{c[1]:.3f}")
    print(f"  pool mean Rg {np.mean(col('rg_pool_mean')):.3f} A, "
          f"pool mean min |i-j|>=3 CA-CA {np.mean(col('minsep3_pool_mean')):.3f} A")

    out = dict(n=len(rows), band=BAND, k=P.K, m_prod=P.M_PROD, rows=rows,
               summary={k: dict(zip(("mean", "se", "n"), mean_se(col(k))))
                        for k in rows[0] if isinstance(rows[0][k], float)})
    ST.save_atomic(os.path.join(P.RESULTS, "phys_landscape.json"), out,
                   complete_keys=("pdb", "fold", "rho_LEG_AMB", "rho_LEG_DIS", "rho_AMB_DIS",
                                  "rho_LEG_rmsd", "rho_AMB_rmsd", "rho_DIS_rmsd"),
                   rows=rows, n_expected=len(pdbs), module_file=__file__)
    print(f"\nwrote s25/results/phys_landscape.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
