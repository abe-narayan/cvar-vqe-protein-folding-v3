"""SPRINT 14, ENER-0 -- what IS the cached AMBER subset, and what is it conditioned on?

The single most valuable asset in the repository carries genuine AMBER on only 2,955 of
262,144 configurations per target.  Before any AMBER statistic is quoted, that subset's
selection has to be established, because if it is not uniform then EVERY AMBER number is
conditioned on the selection -- and one of its three strata is selected by TRUE RMSD.

From `s13/qarch_enum.py` the construction is:
    1,200 uniform indices                                 kind 0
    1,200 indices sampled from the leakage-safe 1-local prior   kind 1
      600 indices drawn from the 2,621 lowest-TRUE-RMSD configurations   kind 2  (ORACLE)
        + the ORACLE snap index
    de-duplicated by np.unique.

This module verifies that claim against the data, quantifies the resulting distortion,
and pins the rule for the rest of the sprint.

    python -m s14.ener_subset
"""
from __future__ import annotations

import numpy as np

from s14 import ener_lib as E


def audit_one(pdb):
    z = E.enum(pdb)
    R = z.rmsd
    sub = R[z.amber_idx]
    rows = {}
    for kind, name in ((0, "uniform"), (1, "prior"), (2, "oracle_band")):
        m = z.amber_kind == kind
        rr = sub[m]
        rows[name] = dict(n=int(m.sum()), rmsd_mean=float(rr.mean()),
                          rmsd_min=float(rr.min()), rmsd_max=float(rr.max()),
                          # where does this stratum sit in the FULL population?
                          pctile_mean=float(np.mean([(R < v).mean() for v in rr[:400]])))
    # the band's defining threshold
    band = sub[z.amber_kind == 2]
    thr = float(band.max())
    n_below = int((R <= thr).sum())
    return dict(
        pdb=pdb, n_subset=int(len(z.amber_idx)), n_full=int(z.B),
        frac=float(len(z.amber_idx) / z.B),
        full_rmsd_mean=float(R.mean()), subset_rmsd_mean=float(sub.mean()),
        subset_minus_full=float(sub.mean() - R.mean()),
        uniform_rmsd_mean=float(sub[z.amber_kind == 0].mean()),
        uniform_minus_full=float(sub[z.amber_kind == 0].mean() - R.mean()),
        band_threshold_rmsd=thr,
        band_pop_frac=float(n_below / z.B),
        strata=rows,
        # is the uniform stratum really uniform in RMSD?  KS-style max deviation
        ks_uniform=float(_ks(sub[z.amber_kind == 0], R)),
        ks_all=float(_ks(sub, R)),
        ks_prior=float(_ks(sub[z.amber_kind == 1], R)),
    )


def _ks(sample, population):
    """Max |F_sample - F_population| on the population's own grid."""
    p = np.sort(population)
    s = np.sort(sample)
    grid = p[np.linspace(0, len(p) - 1, 400).astype(int)]
    Fp = np.searchsorted(p, grid, "right") / len(p)
    Fs = np.searchsorted(s, grid, "right") / len(s)
    return np.abs(Fp - Fs).max()


def main():
    rows = [audit_one(p) for p in E.ENUM_TARGETS]
    for r in rows:
        print(f"{r['pdb']}  subset {r['n_subset']} ({r['frac']*100:.2f}% of space)  "
              f"full_mean {r['full_rmsd_mean']:.3f}  subset_mean {r['subset_rmsd_mean']:.3f} "
              f"({r['subset_minus_full']:+.3f})  uniform_mean {r['uniform_rmsd_mean']:.3f} "
              f"({r['uniform_minus_full']:+.3f})  KS_all {r['ks_all']:.3f}  "
              f"KS_uniform {r['ks_uniform']:.3f}  band_thr {r['band_threshold_rmsd']:.3f} A "
              f"= best {r['band_pop_frac']*100:.2f}%", flush=True)
    agg = dict(
        subset_minus_full=float(np.mean([r["subset_minus_full"] for r in rows])),
        uniform_minus_full=float(np.mean([r["uniform_minus_full"] for r in rows])),
        ks_all=float(np.mean([r["ks_all"] for r in rows])),
        ks_uniform=float(np.mean([r["ks_uniform"] for r in rows])),
        band_pop_frac=float(np.mean([r["band_pop_frac"] for r in rows])),
    )
    print("\nPOOLED:", agg, flush=True)
    E.write("ener_subset", dict(
        what="provenance and distortion of the cached AMBER subset",
        rule=("kind==0 (uniform, ~1194/target) is the ONLY unconditioned AMBER stratum. "
              "kind==2 is selected by TRUE RMSD and is an ORACLE population."),
        pooled=agg, per_target=rows), n_expected=len(E.ENUM_TARGETS))
    return rows


if __name__ == "__main__":
    main()
